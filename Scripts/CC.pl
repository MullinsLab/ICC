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
use lib "/home/wdeng/ICC_v1.1/Scripts/lib";
use seqAlign;
use Getopt::Long;
use File::Path;

my %option = (
	'if' => '',
	'of' => '',
	'cf' => 0.05,
	'mt' => 10,
	'mm' => -10,
	'gp' => -10,
);

my $usage = "\nusage: perl CC.pl [-option value]

options:  
-if		input sequence fasta file
-of		output corrected sequence fasta file
-cf		carry forward cutoff (default: $option{'cf'}, carry forward whose frequency is fewer than cutoff will be corrected)
-mt		match score (default: $option{'mt'})
-mm		mismatch score (default:  $option{'mm'})
-gp		gap penalty (default: $option{'gp'})
		
";

GetOptions (\%option, 'if=s', 'of=s', 'cf=f', 'mt=f', 'mm=f', 'gp=f');

my $inFile = $option{'if'} or die $usage;
my $outFile  = $option{'of'} or die $usage;
my $cfCut = $option{'cf'};
my $match = $option{'mt'};
my $misMatch = $option{'mm'};
my $gapPenalty = $option{'gp'};

my @flowNas = ('T', 'A', 'C', 'G');
my %flowNaHash = (
	'T' => 0,
	'A' => 1,
	'C' => 2,
	'G' => 3,
);
my $name = '';
my (@names, %nameSeq, %picked);
my $readCount = 0;
open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		push @names, $name;
		if ($name =~ /(\d+)$/) {
			$readCount += $1;
		}
	}else {
		#$line =~ s/\-//g;
		$nameSeq{$name} .= uc $line;
	}
}
close IN;

my $idx = my $count = 0;
my $outputDir = "./CC_otus";
if (-e $outputDir) {
	rmtree($outputDir);
}
mkdir $outputDir;
open OUTPUT,">$outFile" or die "couldn't open $outFile: $!\n";
while (@names) {
	my $central_name = shift @names;
	unless ($picked{$central_name}) {	# central_name hasn't been picked by previous clustering
		$idx++;
		$count++;
		my @cluster_names = ();
		my $central_seq = $nameSeq{$central_name};
		print OUTPUT ">$central_name\n$central_seq\n";
		#$central_seq =~ s/\-//g;
		
		push @cluster_names, $central_name;		
		my $duplicates = 1;
		if ($central_name =~ /_(\d+)$/) {
			$duplicates = $1;
		}
		my $readFreq = $duplicates / $readCount;
		if ($readFreq >= 0.01) {	# frequency is above 1%
			foreach my $rest_name (@names) {
				unless ($picked{$rest_name}) {	# rest name hasn't been picked by previous clustering
					my $rest_dup = 1;
					if ($rest_name =~ /_(\d+)$/) {
						$rest_dup = $1;
					}
					my $restFreq = $rest_dup / $readCount;
					my $restPercent = $restFreq / ($readFreq + $restFreq);
					if ($restPercent < $cfCut) {	# carry forward error frequency is as largist as 5% of total frequency in the cluster
						my $rest_seq = $nameSeq{$rest_name};					
						#$rest_seq =~ s/\-//g;
						
						my $etStrs = seqAlign::PairwiseAlign ($central_seq, $rest_seq, $match, $misMatch, $gapPenalty);
						
						if (@$etStrs) {
							foreach my $etStr (@$etStrs) {
								unless ($etStr =~ /DD+/) {
									if ($etStr =~ /D/ && $etStr =~ /I/) {
										#print "$etStr\n";
										my @alignCentralNas = my @alignRestNas = ();
										my $d_char = my $i_char = '';
										my $alignLen = length $etStr;
										my @centralNas = split //, $central_seq;
										my @restNas = split //, $rest_seq;
										my @ets = split //, $etStr;
										my $indelIdx = 0;
										my $insStretch = my $delStretch = my $indels = 0;
										my (%alignIndelIdx, %ins, %del, %indelAlignIdx);
										for (my $i = 0; $i < $alignLen; $i++) {
											if ($ets[$i] eq 'I') {
												++$indelIdx;
												$i_char = shift @restNas;
												$ins{$indelIdx} = $i_char;
												$alignIndelIdx{$i} = $indelIdx;
												$indelAlignIdx{$indelIdx} = $i;
												push @alignCentralNas, '-';
												push @alignRestNas, $i_char;
												if (!$insStretch) {
													$insStretch = 1;
													$indels++;
												}
											}elsif ($ets[$i] eq 'D') {
												++$indelIdx;
												$d_char = shift @centralNas;
												$del{$indelIdx} = $d_char;
												$alignIndelIdx{$i} = $indelIdx;
												$indelAlignIdx{$indelIdx} = $i;
												push @alignCentralNas, $d_char;
												push @alignRestNas, '-';
												if (!$delStretch) {
													$delStretch = 1;
													$indels++;
												}
											}else {
												push @alignCentralNas, shift @centralNas;
												push @alignRestNas, shift @restNas;
												$insStretch = $delStretch = 0;
											}
										}
																			
										my $flag = 0;
										for (my $i = 0; $i < $alignLen; $i++) {
											if ($alignRestNas[$i] eq '-') {	# single deletion
												my $indelIdx = $alignIndelIdx{$i};
												my $d_char = $del{$indelIdx};
												if ($alignCentralNas[$i+1] && $alignCentralNas[$i+2] && $alignCentralNas[$i] eq $alignCentralNas[$i+1] && $alignCentralNas[$i] eq $alignCentralNas[$i+2]) {
													if ($ins{$indelIdx+1}) {
														my $betweenNas = my $pattern = '';
														for (my $j = $i+1; $j < $indelAlignIdx{$indelIdx+1}; $j++) {
															$betweenNas .= $alignRestNas[$j];
														}
														for (my $k = $flowNaHash{$d_char}; $k <= $flowNaHash{$d_char}+4; $k++) {
															my $m = $k % 4;
															$pattern .= $flowNas[$m] . '*';
														}
														if ($betweenNas =~ /^$pattern$/ && $ins{$indelIdx+1} eq $d_char) {
															$flag = 1;
															last;
														}
													}
													if ($ins{$indelIdx-1}) {
														my $betweenNas = my $pattern = '';
														for (my $j = $indelAlignIdx{$indelIdx-1}+1; $j < $i; $j++) {
															$betweenNas .= $alignRestNas[$j];
														}
														for (my $k = $flowNaHash{$d_char}; $k <= $flowNaHash{$d_char}+4; $k++) {
															my $m = $k % 4;
															$pattern .= $flowNas[$m] . '*';
														}
														if ($betweenNas =~ /^$pattern$/ && $ins{$indelIdx-1} eq $d_char) {
															$flag = 1;
															last;
														}
													}
												}elsif ($i > 1 && $alignCentralNas[$i-1] && $alignCentralNas[$i-2] && $alignCentralNas[$i] eq $alignCentralNas[$i-1] && $alignCentralNas[$i] eq $alignCentralNas[$i-2]) {
													if ($ins{$indelIdx+1}) {
														my $betweenNas = my $pattern = '';
														for (my $j = $i+1; $j < $indelAlignIdx{$indelIdx+1}; $j++) {
															$betweenNas .= $alignRestNas[$j];
														}
														for (my $k = $flowNaHash{$d_char}; $k <= $flowNaHash{$d_char}+4; $k++) {
															my $m = $k % 4;
															$pattern .= $flowNas[$m] . '*';
														}
														if ($betweenNas =~ /^$pattern$/ && $ins{$indelIdx+1} eq $d_char) {
															$flag = 1;
															last;
														}
													}
													if ($ins{$indelIdx-1}) {
														my $betweenNas = my $pattern = '';
														for (my $j = $indelAlignIdx{$indelIdx-1}+1; $j < $i; $j++) {
															$betweenNas .= $alignRestNas[$j];
														}
														for (my $k = $flowNaHash{$d_char}; $k <= $flowNaHash{$d_char}+4; $k++) {
															my $m = $k % 4;
															$pattern .= $flowNas[$m] . '*';
														}
														if ($betweenNas =~ /^$pattern$/ && $ins{$indelIdx-1} eq $d_char) {
															$flag = 1;
															last;
														}
													}
												}
											}
										}
										
										if ($flag) {
											$picked{$rest_name} = 1;
											push @cluster_names, $rest_name;
											$count++;
											print OUTPUT ">$rest_name\n$nameSeq{$central_name}\n";
											last;
										}																													
									}
								}							
							}
						}
					}					
				}
			}
		}		
		my $localOutFile = "$outputDir/$idx.fas";
		open OUT, ">$localOutFile" or die "couldn't open $localOutFile: $!\n";
		foreach my $name (@cluster_names) {
			print OUT ">$name\n$nameSeq{$name}\n";
		}
		close OUT;
	}	
}

close OUTPUT;
print "Done!\n";
