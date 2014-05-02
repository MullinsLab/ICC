#!/usr/bin/perl -w

#######################################################################################
# Copyright © 2013 Mullins Lab, University of Washington
# SOFTWARE COPYRIGHT NOTICE 
# This program and its documentation are the copyright of the Mullins Lab, University
# of Washington. All rights are reserved.
# 
# This program is free software. It is supplied without any warranty. You can 
# redistribute it and/or modify it. Mullins Lab is not responsible for its use, 
# misuse, or functionality.
#
# This program was developed by Wenjie Deng <dengw@uw.edu>
#######################################################################################

use strict;
use Getopt::Long;
use File::Basename;

my %option = (
	'is' => '',
	'iq' => '',
	'os' => '',
	'oq' => '',
	'l'  => 100,
	'q'  => 25,
	'h'  => '', 
);

my $usage = "\nUsage: perl readQualFilter.pl [-option value]

options:  
-is     input reads fasta file
-iq     input quality file
-os     output reads fasta file
-oq     output quality file
-l      length cutoff (default: $option{'l'})
-q      average quality score cutoff (default: $option{'q'})
-h      show usage help

";

GetOptions (\%option, 'is=s', 'iq=s', 'os=s', 'oq=s', 'l=i','q=i', 'h');

my $inFasta = $option{'is'} or die $usage;
my $inQual = $option{'iq'} or die $usage;
my $outFasta = $option{'os'} or die $usage;
my $outQual = $option{'oq'} or die $usage;
my $lenCut = $option{'l'};
my $qualCut = $option{'q'};
my $help = $option{'h'};
die $usage if $help;
unless ($lenCut =~ /^\d+$/) {
	die "length cutoff must be a positive integer number\n";
}
unless ($qualCut =~ /^\d+$/) {
	die "average quality cutoff must be a positive integer number\n";
}

my $totalSeq = my $totalQual = 0;
my $seq = my $seqName = "";
my (%nameStatus, %passFilterName);
my $belowLenCutCount = my $belowAllowNCount = my $belowQualCutCount = my $passFilterCount = my $failedCount = 0;
open INFASTA, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <INFASTA>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		if ($seq) {
			my $passFilterFlag = 1;
			my $len = length $seq;
			if ($len < $lenCut) {
				$passFilterFlag = 0;
				$belowLenCutCount++;
			}elsif ($seq =~ /N/) {	# reads with N
				$passFilterFlag = 0;
				$belowAllowNCount++;					
			}
			if ($passFilterFlag) {
				$passFilterName{$seqName} = 1;
			}
		}
		$seqName = $1;
		$seq = '';
		$totalSeq++;
	}else {
		$seq .= uc $line;
	}
}
# last sequence
if ($seq) {
	my $passFilterFlag = 1;
	my $len = length $seq;
	if ($len < $lenCut) {
		$passFilterFlag = 0;
		$belowLenCutCount++;
	}elsif ($seq =~ /N/) {	# reads with N		
		$passFilterFlag = 0;
		$belowAllowNCount++;
	}
	if ($passFilterFlag) {
		$passFilterName{$seqName} = 1;
	}
}
$seqName = $seq = '';
close INFASTA;

print "total $totalSeq in seq file.\n$belowLenCutCount reads length shorter than $lenCut.\n";
print "$belowAllowNCount reads contain ambiguities.\n";

if ($qualCut > 0) {
	my $belowQualCutCount = 0;
	my @quals = ();
	open INQUAL, $inQual or die "couldn't open $inQual: $!\n";
	while (my $line = <INQUAL>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			if (@quals) {
				my $totalQuals = 0;
				foreach my $qual (@quals) {
					$totalQuals += $qual;
				}
				my $avgQual = $totalQuals / scalar @quals;
				if ($avgQual < $qualCut) {
					$passFilterName{$seqName} = 0;
					$belowQualCutCount++;
				}
			}
			$seqName = $1;
			@quals = ();
			$totalQual++;
		}elsif ($passFilterName{$seqName}) {
			my @partialQuals = split /\s+/, $line;
			push @quals, @partialQuals;
		}
	}
	if (@quals) {
		my $totalQuals = 0;
		foreach my $qual (@quals) {
			$totalQuals += $qual;
		}
		my $avgQual = $totalQuals / scalar @quals;
		if ($avgQual < $qualCut) {
			$passFilterName{$seqName} = 0;
			$belowQualCutCount++;
		}
	}
	$seqName = '';
	close INQUAL;
	print "total $totalQual in qual file\n";
	print "$belowQualCutCount reads average quality score below $qualCut.\n";
}

open OUTFASTA, ">$outFasta" or die "couldn't open $outFasta: $!\n";
open INFASTA, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <INFASTA>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$seqName = $1;
		if ($passFilterName{$seqName}) {
			print OUTFASTA ">$seqName\n";
			$passFilterCount++;
		}else {
			$failedCount++;
		}
	}elsif ($passFilterName{$seqName}) {
		print OUTFASTA "$line\n";
	}
}
close OUTFASTA;
close INFASTA;
$seqName = '';

open OUTQUAL, ">$outQual" or die "couldn't open $outQual: $!\n";
open INQUAL, $inQual or die "couldn't open $inQual: $!\n";
while (my $line = <INQUAL>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$seqName = $1;
		if ($passFilterName{$seqName}) {
			print OUTQUAL ">$seqName\n";
		}
	}elsif ($passFilterName{$seqName}) {
		print OUTQUAL "$line\n";
	}
}
close OUTQUAL;
close INQUAL;

print "combined there are $passFilterCount meet the cutoff. $failedCount don't.\n";
print "All done.\n";

