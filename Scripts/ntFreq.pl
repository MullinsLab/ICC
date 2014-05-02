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

my $usage = "\nusage: perl ntFreq.pl [-option value]
options:  
-ia		input alignment file
-of		output frequency file
-uf		unique sequence flag (default: false). indicating if the sequences in alignment are unique,
		if true, the duplicate info must be at the end of sequence name
		
";

my $inFile = '';
my $FreqFile = '';
my $uniqFlag = '';

GetOptions ('ia=s' => \$inFile, 'of=s' => \$FreqFile, 'uf' => \$uniqFlag);
die $usage unless ($inFile && $FreqFile);

my $flag = my $totalSeq = my $refFlag = 0;
my $duplicates = 1;
my $refName = my $refSeq = '';
my ($seqName, $seq, $len, $nameSeq, $naStart, $naEnd, @seqNames, @fwdseqNames, @revseqNames, %seqDuplicate);
die "The input file is empty.\n" if -z $inFile;
open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	if ($line =~ />Reference/i) {
		$refName = "Reference_1";	# indicates one duplicate for the reference
		$refFlag = 1;
	}elsif ($refFlag) {
		if ($line =~ /^>/) {
			$refFlag = 0;
		}else {
			$refSeq .= $line;
		}		
	}
	if (!$refFlag) {
		if ($line =~ /^>(.*)$/) {
			if ($flag) {
				my $seqLen = length $seq;
				if ($flag == 1) {
					$len = $seqLen;
				}else {
					if ($seqLen != $len) {
						die "error: length of sequences are not same, your sequences are probably not aligned.\n";
					}
				}
				my @seqNas = split //, $seq;
				$nameSeq->{$seqName} = \@seqNas;
			}
			$seqName = $1;
			push @seqNames, $seqName;
			$seq = "";
			$flag++;
			if ($uniqFlag) {	# there is duplicates info at the end of sequence name
				if ($seqName =~ /(\d+)$/) {
					$duplicates = $1;
				}else {
					die "there is no duplicates info at the end of sequence name\n";
				}
			}else {
				$duplicates = 1;
			}			
			$seqDuplicate{$seqName} = $duplicates;
			$totalSeq += $duplicates;
		}else {
			$seq .= uc $line;
		}
	}
}
if ($flag) {
	my @seqNas = split //, $seq;
	my $seqLen = length $seq;
	if ($flag > 1 && $seqLen != $len) {
		die "error: length of sequences are not equal, your sequences are probably not aligned.\n";
	}
	$nameSeq->{$seqName} = \@seqNas;
}
close IN;

open STAT, ">$FreqFile" or die "couldn't open $FreqFile: $!\n";
#print STAT "There are total $totalSeq sequences. alignment length $len.\n\n";
print STAT "Position\tConsensus";

if ($nameSeq) {
	my $posNaCount = CountFreq (\@seqNames, $nameSeq, \%seqDuplicate, $len);
	WriteFreq ($totalSeq, $posNaCount, $len);
}

close STAT;
print "Done!\n";


sub CountFreq {
	my ($sub_seqNames, $sub_nameSeq, $sub_seqDuplicate, $sub_len) = @_;
	my ($posNaCount);
	for (my $i = 0; $i < $sub_len; $i++) {
		foreach my $seqName (@$sub_seqNames) {			
			if (!$posNaCount->{$i}->{$sub_nameSeq->{$seqName}->[$i]}) {
				$posNaCount->{$i}->{$sub_nameSeq->{$seqName}->[$i]} = 0;
			}
			$posNaCount->{$i}->{$sub_nameSeq->{$seqName}->[$i]} += $sub_seqDuplicate->{$seqName};											
		}
	}
	return $posNaCount;
}

sub WriteFreq {
	my ($sub_totalSeq, $sub_posNaCount, $sub_len) = @_;
	my @nts = ('-', 'A', 'C', 'G', 'T');
	foreach my $nt (@nts) {
		print STAT "\t$nt";
	}
	print STAT "\tReads\n";
	for (my $i = 0; $i < $sub_len; $i++) {
		my $pos = $i + 1;
		print STAT "$pos";
		my $flag = 0;
		my $cons = '';
		foreach my $nt (@nts) {
			$sub_posNaCount->{$i}->{$nt} = 0 if (!$sub_posNaCount->{$i}->{$nt});
			if (!$flag) {
				$cons = $nt;
				++$flag;
			}else {
				if ($sub_posNaCount->{$i}->{$nt} >= $sub_posNaCount->{$i}->{$cons}) {
					$cons = $nt;
				}
			}
		}
		print STAT "\t$cons";
		foreach my $nt (@nts) {
			if ($sub_totalSeq) {
				if ($sub_posNaCount->{$i}->{$nt}) {
					my $percent = int ($sub_posNaCount->{$i}->{$nt} / $sub_totalSeq * 1000000 + 0.5) / 1000000;				
					print STAT "\t$percent";				
				}else {
					print STAT "\t0";
				}
			}else {
				print STAT "\t0";
			}			
		}
		if ($sub_totalSeq) {
			print STAT "\t", $sub_totalSeq;
		}else {
			print STAT "\t0";
		}
		print STAT "\n";
	}
	print STAT "\n";
}
