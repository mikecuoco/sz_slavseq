#!/usr/bin/env perl
$^W=1;

use strict;
use Bio::DB::HTS;
use Data::Dumper;
use Cwd "abs_path";


sub uniqid
{
    return time()."_".int(rand(1e6))
}


sub use_tmp_dir
{
    my $tmpbase=shift;
    #my $tmpbase='/tmp';
    my $old=`pwd`;
    chomp $old;
    my $tmp= "$tmpbase/tmp_".uniqid();
    mkdir $tmp or die "$! \"$tmp\"";
    chdir $tmp or die "$! \"$tmp\"";
    $tmp=`pwd`;
    chomp $tmp;
    return ($old, $tmp);
}


sub remove_tmp_and_return_to_old_dir
{
    my ($old,$tmp)=@_;
    die "safety net on tmp dir" unless $tmp =~ /tmp_[a-zA-Z0-9-_]+$/;
    chdir $old or die;
    system("\\rm -fr $tmp")==0 or die;
}


sub prio_pair_rmdup
{
    my ($fn, $ofn)=@_; # store arguments as variables
    my $bam = Bio::DB::HTSfile->open($fn); # read input bam file
  
    #my $good_bam = Bio::DB::HTSfile->open($good_ofn, 'w');
    #my $bad_bam = Bio::DB::HTSfile->open($bad_ofn, 'w');

    #assume reads are sorted by name
    # my $header = $bam->header(); # Bio:DB:Sam
    my $header = $bam->header_read; # read header info
 
    #$good_bam->header_write($header);
    #$bad_bam->header_write($header);
    
    my $target_names = $header->target_name; # return a reference to an array of reference sequence names

	# write to second argument
    open OUT, " | sort -S 5000M -k 4,4 -k 2,2rn -k 3,3rn | uniq -f 3 -c | perl -pe 's/^ +(\\d+) +(\\S+)/\$2\\tXD:i:\$1/' | cut -f 1,2 | sort -S 5000M > $ofn";

    #open OUT, " > $ofn";

    my $r1=undef;
    my $r2=undef;
    
    while($r2 = $bam->read1($header)) # read one alignment from file and return Alignment object
    {
	# filter out secondary matches (bwa mem compatibility)
	# secondary matches won't be used in selection
	#	256 = not primary alignment
	#	2048 = supplementary alignment
	# see https://broadinstitute.github.io/picard/explain-flags.html for SAM flags info
	next if ($r2->flag & 256) || ($r2->flag & 2048); # continue to next iteration if true; flag is SAM bitwise flag field
	

	if(defined($r1))
	{
	    if(
		($r2->qname ne $r1->qname) || # if read names are not the same...
		(($r1->flag & 4) && ($r1->flag & 64)) # drop if R1 unmapped. 4 = unmapped, 64 = first in pair
		)
		
	    {
		#$bad_bam->write1($r1);
	    }
	    else
	    {
		#we got a good pair here.

		if($r1->flag & 128) # 128 = second in pair
		{
			# switch values of r1 and r2
		    my $tmp=$r1;
		    $r1=$r2;
		    $r2=$tmp;
		}

		# quality scores are at each base, so summing across bases of the read
		my $sumqual=0;
		for my $q ($r1->qscore, $r2->qscore)
		{
		    $sumqual+=$q; # sum  of quality scores
		}
		
		# pos is 0-based leftmost coordinate of the aligned sequence on the reference sequence, calend is rightmost
		# 16 = read reverse strand
		my $r1pos=($r1->flag & 16)?$r1->calend+1:$r1->pos+1; # convert from 0-based to 1-based coordinates
		my $r1strand=($r1->flag & 16)?'-':'+'; # marking strand
		
		#my $r2pos=($r2->flag & 4)?'*':($r2->flag & 16)?$r2->calend+1:$r2->pos+1;
		#my $r2strand=($r2->flag & 4)?'*':($r2->flag & 16)?'-':'+';
		
		# join target ID, position, and strand info
		my $pos=join(':',$target_names->[$r1->tid], $r1pos, $r1strand);
		
		# join and print more metrics to output file
		print OUT join("\t", $r1->qname, $r1->qual + $r2->qual, $sumqual, $pos),"\n";
		
		$r2=undef;
	    }
	}

	$r1=$r2;
    }
    close OUT;
}

my $input_bam_fn=shift; 
my $output_bam_fn=shift;

my $input_bam_fn_abspath=abs_path($input_bam_fn);
die "Output file exists: $output_bam_fn" if -e $output_bam_fn;

$ENV{LC_ALL}="C";
my $pwd=`pwd`;
chomp $pwd;
my ($curdir, $tmpdir)=use_tmp_dir $pwd;
chdir $curdir;
system("ln -s $input_bam_fn_abspath $tmpdir/input.bam");
chdir $tmpdir;

prio_pair_rmdup 'input.bam', 'selected.txt';

system("samtools view -H input.bam > header.txt"); # get header only
# sort input bam by first column, output lines identical to selected.txt, and append header to create output bam file 
system(qq[samtools view input.bam | sort -T ./ -S 1500M -s -k 1,1 | join -t "\t" - selected.txt | cat header.txt - | samtools view -S -b - > output.bam]);
chdir $curdir;
rename "$tmpdir/output.bam", "$output_bam_fn";

remove_tmp_and_return_to_old_dir($curdir, $tmpdir);
