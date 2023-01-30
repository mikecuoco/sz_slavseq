#!/usr/bin/env python
__author__ = "Apuã Paquola"


def overlaps(iv1, iv2):
    """Returns a boolean indicating whether intervals overlap."""
    return (
        iv1.chrom == iv2.chrom and iv1.end - 1 >= iv2.start and iv2.end - 1 >= iv1.start
    )


def is_within(iv1, iv2):
    """Returns a boolean indicating whether iv1 is within iv2."""
    return iv1.chrom == iv2.chrom and iv1.end <= iv2.end and iv1.start >= iv2.start


def intersection(iv1, iv2):
    """Returns the overlapping part of two intervals if they overlap or
    None if they don't.
    """
    if overlaps(iv1, iv2):
        return Interval(iv1.chrom, max(iv1.start, iv2.start), min(iv1.end, iv2.end))
    else:
        return None


class Interval:
    """chrom is a string, start is 0-based, end is 1-based"""

    def __init__(self, chrom=None, start=None, end=None):
        assert start <= end
        self.chrom = chrom
        self.start = start
        self.end = end

    def as_tuple(self):
        return self.chrom, self.start, self.end

    def __str__(self):
        return str(self.chrom) + ":" + str(self.start) + "-" + str(self.end)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __lt__(self, other):
        return self.as_tuple() < other.as_tuple()

    def __hash__(self):
        return self.as_tuple().__hash__()
