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

my %option = (
	'if' => '',
	'of' => '',
	'nt' => 'UniqRead_',
	'uf' => '',
	'df' => ''
);

my $usage = "\nusage: perl uniqueReads.pl [-option value]

options:  
-if		input fasta file
-of		output fasta file
-nt		name tag (default: 'uniqRead_')
-uf		unique flag (default: false). Identify if the input file has duplicate info at the end of sequence name
-df		direction flag (default: false). If set, find unique sequences based on forward (names start with F_) 
		and reverse (names start with R_), if not, find unique sequences based on whole data set

";

GetOptions (\%option, 'if=s', 'of=s', 'nt=s', 'sf=s', 'uf', 'df');

my $inFile = $option{'if'} or die $usage;
my $outFile = $option{'of'} or die $usage;
my $nameFile = $outFile.".name";
my $nameTag = $option{'nt'};
my $uniqFlag = $option{'uf'};
my $dirFlag = $option{'df'};

my $count = my $fwdCount = my $revCount = my $refFlag = 0;
my $seq = my $seqName = my $refName = my $refSeq = '';
my $duplicates = 1;
my (%seqCount, %fwdSeqCount, %revSeqCount);
my %seqLen;
my %nameTag;
my (@uniqSeqs, @fwdUniqSeqs, @revUniqSeqs, $uniqDup, $fwdUniqDup, $revUniqDup);
open INFASTA, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <INFASTA>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	if ($line =~ /Reference/i) {
		$refName = "Reference_1";	# indicates one duplicate for the reference
		$refFlag = 1;
	}elsif ($refFlag) {
		if ($line =~ /^>/) {
			$refFlag = 0;
		}else {
			$refSeq .= $line;
			print "ref: $refSeq\n";
		}		
	}
	if (!$refFlag) {
		if ($line =~ /^>(\S+)/) {
			if ($seq) {
				unless ($seq =~ /^\-+$/) {	# in case all gaps in sequence
					if ($dirFlag) {
						if ($seqName =~ /^F_/) {
							if (!$fwdSeqCount{$seq}) {
								$fwdSeqCount{$seq} = 0;
								push @fwdUniqSeqs, $seq;
							}
							push @{$fwdUniqDup->{$seq}}, $seqName;
							$fwdSeqCount{$seq} += $duplicates;	
							$fwdCount += $duplicates;
						}elsif ($seqName =~ /^R_/) {
							if (!$revSeqCount{$seq}) {
								$revSeqCount{$seq} = 0;
								push @revUniqSeqs, $seq;
							}
							push @{$revUniqDup->{$seq}}, $seqName;
							$revSeqCount{$seq} += $duplicates;	
							$revCount += $duplicates;
						}else {
							die "No direction information in $seqName\n";
						}
					}else {
						if (!$seqCount{$seq}) {
							$seqCount{$seq} = 0;
							push @uniqSeqs, $seq;
						}
						push @{$uniqDup->{$seq}}, $seqName;
						$seqCount{$seq} += $duplicates;							
					}
					$count += $duplicates;					
				}		
				$seqName = $seq = "";
			}
			$seqName = $1;
			if ($uniqFlag) {	# there is duplicates info at the end of sequence name
				if ($seqName =~ /(\d+)$/) {
					$duplicates = $1;
				}else {
					die "there is no duplicates info at the end of sequence name\n";
				}
			}else {
				$duplicates = 1;
			}			
		}else {
			$seq .= uc $line;
		}
	}	
}
if ($seq) {
	unless ($seq =~ /^\-+$/) {	# in case all gaps in sequence
		if ($dirFlag) {
			if ($seqName =~ /^F_/) {
				if (!$fwdSeqCount{$seq}) {
					$fwdSeqCount{$seq} = 0;
					push @fwdUniqSeqs, $seq;
				}
				push @{$fwdUniqDup->{$seq}}, $seqName;
				$fwdSeqCount{$seq} += $duplicates;	
				$fwdCount += $duplicates;
			}elsif ($seqName =~ /^R_/) {
				if (!$revSeqCount{$seq}) {
					$revSeqCount{$seq} = 0;
					push @revUniqSeqs, $seq;
				}
				push @{$revUniqDup->{$seq}}, $seqName;
				$revSeqCount{$seq} += $duplicates;	
				$revCount += $duplicates;
			}else {
				die "No direction information in $seqName\n";
			}
		}else {
			if (!$seqCount{$seq}) {
				$seqCount{$seq} = 0;
				push @uniqSeqs, $seq;
			}
			push @{$uniqDup->{$seq}}, $seqName;
			$seqCount{$seq} += $duplicates;			
		}
		$count += $duplicates;	
	}
	$seqName = $seq = "";
}
close INFASTA;
print "total $count sequences. ";

my $uniqueCount = 0;
open OUT, ">$outFile" or die "couldn't open $outFile: $!\n";
open NAME, ">$nameFile" or die "couldn't open $nameFile: $!\n";
if ($refName && $refSeq) {
	print OUT ">$refName\n$refSeq\n";
}

if ($dirFlag) {
	my $fwdUniqCount = my $revUniqCount = 0;
	foreach my $seq (sort {$fwdSeqCount{$b} <=> $fwdSeqCount{$a}} keys %fwdSeqCount) {
		$fwdUniqCount++;
		my $name = "F_".$nameTag.$fwdUniqCount."_".$fwdSeqCount{$seq};
		print OUT ">",$name,"\n",$seq,"\n";
		print NAME $name, " ", join(',', @{$fwdUniqDup->{$seq}}), "\n";
	}
	foreach my $seq (sort {$revSeqCount{$b} <=> $revSeqCount{$a}} keys %revSeqCount) {
		$revUniqCount++;
		my $name = "R_".$nameTag.$revUniqCount."_".$revSeqCount{$seq};
		print OUT ">",$name,"\n",$seq,"\n";
		print NAME $name, " ", join(',', @{$revUniqDup->{$seq}}), "\n";
	}
	print "$fwdCount forward reads, $fwdUniqCount unique forward reads, $revCount reverese reads, $revUniqCount unique reverse reads\n";
}else {
	foreach my $seq (sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		$uniqueCount++;
		my $name = $nameTag.$uniqueCount."_".$seqCount{$seq};
		print OUT ">",$name,"\n",$seq,"\n";
		print NAME $name, " ", join(',', @{$uniqDup->{$seq}}), "\n";
	}
	print "$uniqueCount unique reads\n";
}
close OUT;
