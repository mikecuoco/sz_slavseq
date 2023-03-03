import pysam
import pandas as pd
import numpy as np
from scipy.stats import poisson
from collections import deque, namedtuple


class BasePeakCaller(object):
	def __init__(self, bam: pysam.AlignmentFile, read_filter: bool = True) -> None:
		self.bam = bam
		self.read_filter = read_filter

	def get_reads(self, contig):
		"""
		Generator that fetches reads from bam file.
		input:
		:param contig: str, contig to fetch reads from
		output:
		:yeilds: pysam.AlignedSegment
		"""
		Read = namedtuple("Read", ["chr","start","end","is_reverse"])
		for r in self.bam.fetch(contig, multiple_iterators=True):
			if self.read_filter:
				if (
					r.is_read1
					and (not r.is_unmapped)
					and (not r.is_secondary)
					and (not r.is_supplementary)
					# and (not r.has_tag("XA"))
					# and (not r.has_tag("SA"))
					and (r.has_tag("YA") and r.has_tag("YG"))
					and (r.get_tag("YA") > 20 and r.get_tag("YA") > r.get_tag("YG"))
					# and (r.mapping_quality >= 60)
				):
					yield Read(r.reference_name,r.reference_start,r.reference_end,r.is_reverse)
			else:
				yield Read(r.reference_name,r.reference_start,r.reference_end,r.is_reverse)

	def run_peak_caller(self, contigs: list = None):
		"""
		Runs peak caller on all contigs in bam file.
		"""

		if contigs is None:
			contigs = self.bam.references
		else:
			for c in contigs:
				assert c in self.bam.references, f"{c} not in bam file."

		# run peak caller on each contig
		result = [p for c in contigs for p in self.call_peaks(c)]

		return pd.DataFrame(result)


class SlidingPeakCaller(BasePeakCaller):
	def __init__(
		self,
		bam: pysam.AlignmentFile,
		read_filter: bool,
		window_size: int,
		step_size: int,
		min_reads: int,
		merge: bool,
		bg_sizes: list,
	) -> None:
		"""
		Initializes a PeakCaller object with a bam file and peak and background window sizes.
		input:
		:param bam: pysam.AlignmentFile, opened bam file
		:param window_size: int, size of window to find candidate peaks
		:param step_size: int, size of step between windows for candidate peaks
		:param bg_sizes: list, sizes of background windows
		NOTE: Slow with high numbers of bg_sizes, since each size requires bam file to be opened.
		"""
		super().__init__(bam, read_filter)
		if bg_sizes:
			assert window_size < min(
				bg_sizes
			), "peak size must be smaller than background window sizes"
		self.window_size = window_size
		self.step_size = step_size
		self.min_reads = min_reads
		self.merge = merge
		self.bg_sizes = bg_sizes

	def window(
		self, contig: str, window_size: int, step_size: int = 1, min_reads: int = 1
	):
		"""
		Generator that slides window of size window_size across a contig and yields windows with > 0 read alignments.
		:param contig: str, contig in bam
		:param window_size: int, size of window
		:param step_size: int, size of step bewteen windows
		:param min_reads: int, minimum number of reads in window
		"""

		try:
			# get reads for this contig
			f = self.get_reads(contig)
			r = next(f)
		except StopIteration:
			return

		# initialize deque of reads
		# data structure for fast removal of reads outside of sliding window and addition of reads inside of sliding window
		freads, rreads = deque(), deque()

		reflen = self.bam.get_reference_length(contig)
		for start in range(0, reflen - window_size + 1, step_size):
			end = start + window_size if start + window_size < reflen else reflen
			# check if reads need to be removed
			# reads is False if reads is empty
			while freads and freads[0].start < start:
				freads.popleft()
			while rreads and rreads[0].start < start:
				rreads.popleft()

			# check if reads need to be added
			while (r.start < end) and (r.start >= start):
				if r.is_reverse:
					rreads.append(r)
				else:
					freads.append(r)
				try:
					r = next(f)
				except StopIteration:
					if len(freads) + len(rreads) > min_reads:
						return {
							"chr": contig,
							"start": start,
							"end": end,
							"width": end - start,
							"center": int(start + window_size / 2),
							"count_fwd": len(freads),
							"count_rev": len(rreads),
							"reads_fwd": set(freads), 
							"reads_rev": set(rreads),
							"is_reverse": True if len(rreads) > len(freads) else False
						}
					else:
						return

			# yield the window
			if len(freads) + len(rreads) > min_reads:
				yield {
					"chr": contig,
					"start": start,
					"end": end,
					"width": end - start,
					"center": int(start + window_size / 2),
					"count_fwd": len(freads),
					"count_rev": len(rreads),
					"reads_fwd": set(freads), 
					"reads_rev": set(rreads),
					"is_reverse": True if len(rreads) > len(freads) else False
				}

	def find_peaks(self, peak_windows, bg_windows: dict, pval_cutoff=10e-5):
		"""
		Generator that yields peaks from read alignments to contig using poisson test against local background.
		:param peak_windows: iterator of windows of reads from bam file
		:param bg_windows: dict of generators, windows of reads from bam file keyed by window size
		"""
		# initialize first background windows
		bg = {}
		for ws, w in bg_windows.items():
			try:
				bg[ws] = next(w)
			except StopIteration:
				break

		# if no background windows, skip contig
		if not bg:
			return

		while True:
			try:
				# get next peak window
				p = next(peak_windows)
			except StopIteration:
				return

			# adjust window centers over peak center
			# for each window_size, move to next window if: peak center is greater than window center
			# test peak against window with largest lambda
			k = "rev" if p["is_reverse"] else "fwd"
			lambdas = []
			for ws in bg.keys():
				while p["center"] > bg[ws]["center"]:
					try:
						bg[ws] = next(bg_windows[ws])
					except StopIteration:
						break

				if ws > p["width"]:
					assert (
						p["width"] < bg[ws]["width"]
					), "peak {}:{}-{} width must be less than background window {}:{}-{} width".format(
									p["chr"],
						p["start"],
						p["end"],
						bg[ws]["chr"],
						bg[ws]["start"],
						bg[ws]["end"],
					)
					assert (
						p["start"] >= bg[ws]["start"]
					), "peak {}:{}-{} start must be greater than or equal to background window {}:{}-{} start".format(
						p["chr"],
						p["start"],
						p["end"],
						bg[ws]["chr"],
						bg[ws]["start"],
						bg[ws]["end"],
					)
					assert (
						p["end"] <= bg[ws]["end"]
					), "peak {}:{}-{} end must be less than or equal to background window {}:{}-{} end".format(
						p["chr"],
						p["start"],
						p["end"],
						bg[ws]["chr"],
						bg[ws]["start"],
						bg[ws]["end"],
					)
					lambdas.append(int(round((bg[ws][f"count_{k}"] / (bg[ws]["end"] - bg[ws]["start"])) * p["width"])))

			# calculate p-value
			mu = max(lambdas)
			p["p"] = poisson._pmf(p[f"count_{k}"], mu)
			p["fold_change"] = (
				np.float16(p[f"count_{k}"] / mu) if mu > 0 else np.float16(p[f"count_{k}"])
			)

			# yield peak
			if p["p"] < pval_cutoff:
				yield p

	def merge_peaks(self, peaks):
		"""
		Merge overlapping peaks
		:param peaks: generator, peaks from bam file, identified by find_peaks
		"""
		# merge the peaks
		try:
			p = next(peaks)  # grab first peak
		except StopIteration:  # skip if no peaks
			return

		while True:
			# if next peak overlaps with current peak, extend current peak,
			# else, yield current peak and set next peak as current peak
			try:
				n = next(peaks)
				if (n["start"] >= p["start"] and n["start"] <= p["end"]) or (
					n["end"] >= p["start"] and n["end"] <= p["end"]
				):
					p["end"] = n["end"]
					p["width"] = p["end"] - p["start"]
					p["center"] = int(p["start"] + p["width"] / 2)
					p["reads_fwd"].union(n["reads_fwd"])
					p["reads_rev"].union(n["reads_rev"])
				else:
					p["count_fwd"] = int(len(p["reads_fwd"]))
					p["count_rev"] = int(len(p["reads_rev"]))
					p["is_reverse"] = True if p["count_rev"] > p["count_fwd"] else False
					yield p
					p = n
			# if no more peaks, return current peak
			except StopIteration:
				p["count_fwd"] = int(len(p["reads_fwd"]))
				p["count_rev"] = int(len(p["reads_rev"]))
				p["is_reverse"] = True if p["count_rev"] > p["count_fwd"] else False
				return p
			

	def call_peaks(self, contig: str):
		"""
		Runs call peaks for a single contig
		:param contig: str, contig to call peaks on
		"""

		# skip small contigs
		if self.window_size > self.bam.get_reference_length(contig):
			return

		# get peak windows
		peak_windows = self.window(contig, self.window_size, self.step_size, min_reads=self.min_reads)

		# adjust for background if bg_sizes is not empty
		if self.bg_sizes:
			bg_windows1 = {ws: self.window(contig, ws) for ws in self.bg_sizes}

			# get significant peaks for largest local background window size
			peaks = self.find_peaks(peak_windows, bg_windows1)

			# merge overlapping peaks
			if self.merge:
				merged = self.merge_peaks(peaks)

				bg_windows2 = {ws: self.window(contig, ws) for ws in self.bg_sizes}

				return self.find_peaks(merged, bg_windows2)
			else:
				return peaks
		else:
			if self.merge:
				return self.merge_peaks(peak_windows)
			else:
				return peak_windows


class OverlapPeakCaller(BasePeakCaller):
	def __init__(self, bam: pysam.AlignmentFile, read_filter: bool) -> None:
		"""
		Initializes a PeakCaller object with a bam file and peak and background window sizes.
		input:
		:param bam: pysam.AlignmentFile, opened bam file
		"""
		super().__init__(bam, read_filter)

	def merge_reads(self, reads):
		"""
		Merge overlapping peaks
		:param peaks: generator, peaks from bam file, identified by find_peaks
		"""
		# merge the peaks
		try:
			p = dict(next(reads)._asdict())  # grab first read
			p["count_fwd"] = 0 if p["is_reverse"] else 1
			p["count_rev"] = 1 if p["is_reverse"] else 0
		except StopIteration:  # skip if no peaks
			return

		while True:
			# if next peak overlaps with current peak, extend current peak,
			# else, yield current peak and set next peak as current peak
			try:
				n = dict(next(reads)._asdict())
				n["count_fwd"] = 0 if n["is_reverse"] else 1
				n["count_rev"] = 1 if n["is_reverse"] else 0
				if (n["start"] >= p["start"] and n["start"] <= p["end"]) or (
					n["end"] >= p["start"] and n["end"] <= p["end"]
				):
					p["end"] = n["end"]
					p["width"] = p["end"] - p["start"]
					p["center"] = int(p["start"] + p["width"] / 2)
					p["count_fwd"] += 0 if n["is_reverse"] else 1
					p["count_rev"] += 1 if n["is_reverse"] else 0
				else:
					p["is_reverse"] = True if p["count_rev"] > p["count_fwd"] else False
					yield p
					p = n
			# if no more reads, return current peak
			except StopIteration:
				p["is_reverse"] = True if p["count_rev"] > p["count_fwd"] else False
				return p
			

	def call_peaks(self, contig: str):

		reads = self.get_reads(contig)
		return self.merge_reads(reads)
