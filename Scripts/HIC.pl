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
use lib "/home/wdeng/ICC_v1.0/Scripts/lib";
use seqAlign;
use utils;
use paths;
use Getopt::Long;
use File::Path;

my %option = (
	'if' => '',
	'of' => '',
	'od' => '',
	'mt' => 10,
	'mm' => -9,
	'gp' => -15,
	'cs' => 2,
);

my $usage = "\nusage: perl HIC.pl [-option value]

options:  
-if	input sequence fasta file
-of	output corrected sequence fasta file
-od	output distance file
-mt	match score (default: $option{'mt'})
-mm	mismatch score (default:  $option{'mm'})
-gp	gap penalty (default: $option{'gp'})
-cs	minimal cluster size as a cluster seed (defualt: $option{'cs'})
		
";

GetOptions (\%option, 'if=s', 'of=s', 'od=s', 'mt=f', 'mm=f', 'gp=f', 'cs=i');

my $inFile = $option{'if'} or die $usage;
my $outFile  = $option{'of'} or die $usage;
my $outDist = $option{'od'} or die $usage;
my $match = $option{'mt'};
my $misMatch = $option{'mm'};
my $gapPenalty = $option{'gp'};
my $minClusterSize = $option{'cs'};
my $scriptPath = $paths::scriptPath;
my $name = '';
my (@names, %nameSeq, %picked);
open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		push @names, $name;
	}else {
		#$line =~ s/\-//g;
		$nameSeq{$name} .= uc $line;
	}
}
close IN;

my $idx = my $count = 0;
my $outputDir = "./HIC_otus";
if (-e $outputDir) {
	rmtree($outputDir);
}
mkdir $outputDir;
open OUTPUT,">$outFile" or die "couldn't open $outFile: $!\n";
open DIST, ">$outDist" or die "couldn't open $outDist: $!\n";
while (@names) {
	my $central_name = shift @names;
	next if ($central_name =~ /^Reference_1$/i);
	unless ($picked{$central_name}) {	# central_name hasn't been picked by previous clustering
		$idx++;
		$count++;
		my @cluster_names = ();
		my $central_seq = $nameSeq{$central_name};
		$central_seq =~ s/\-//g;
		
		push @cluster_names, $central_name;		
		my $central_duplicates = 1;
		if ($central_name =~ /_(\d+)$/) {
			$central_duplicates = $1;
		}
		my $total_rest = 0;
		if ($central_duplicates >= $minClusterSize) {
			foreach my $rest_name (@names) {
				unless ($picked{$rest_name}) {	# rest name hasn't been picked by previous clustering
					my $rest_seq = $nameSeq{$rest_name};
					$rest_seq =~ s/\-//g;
					my $rest_duplicates = 1;
					if ($rest_name =~ /_(\d+)$/) {
						$rest_duplicates = $1;
					}
					my $etStrs = seqAlign::PairwiseAlign ($central_seq, $rest_seq, $match, $misMatch, $gapPenalty);
				
					if (@$etStrs) {
						my $flag = my $indels = 0;
						foreach my $etStr (@$etStrs) {
							my @ets = split //, $etStr;
							my @alignCentralNas = my @alignRestNas = ();
							my $alignLen = length $etStr;
							my @centralNas = split //, $central_seq;
							my @restNas = split //, $rest_seq;
							my $insStretch = my $delStretch = $indels = 0;
							foreach my $et (@ets) {
								if ($et eq 'I') {
									push @alignCentralNas, '-';
									push @alignRestNas, shift @restNas;
									if (!$insStretch) {
										$insStretch = 1;
										$indels++;
									}
								}elsif ($et eq 'D') {
									push @alignCentralNas, shift @centralNas;
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
							
							$flag = 0;
							for (my $i = 0; $i < $alignLen; $i++) {
								if ($alignCentralNas[$i] eq '-' && $alignCentralNas[$i+1] && $alignCentralNas[$i+1] eq '-') {	# two insertions stretch
									if ( $alignRestNas[$i] eq $alignRestNas[$i+1]) {	# same insertions
										if ($i == 0) {
											if ($alignRestNas[$i] eq $alignRestNas[$i+2]) {
												$i++; # homopolymer stretch, set next position after the insertion stretch
											}else {
												$flag = 1;
												last;
											}											
										}elsif (($alignRestNas[$i-1] && $alignRestNas[$i] eq $alignRestNas[$i-1]) || ($alignRestNas[$i+2] && $alignRestNas[$i] eq $alignRestNas[$i+2])) {
											$i++; # homopolymer stretch, set next position after the insertion stretch
										}else {
											$flag = 1;
											last;
										}
									}else {
										$flag = 1;
										last;
									}
								}elsif ($alignRestNas[$i] eq '-' && $alignRestNas[$i+1] && $alignRestNas[$i+1] eq '-') {	# two deletions stretch
									if ( $alignCentralNas[$i] eq $alignCentralNas[$i+1]) {	# same deletions
										if ($i == 0) {
											if ($alignCentralNas[$i] eq $alignCentralNas[$i+2]) {
												$i++; # homopolymer stretch, set next position after the deletion stretch
											}else {
												$flag = 1;
												last;
											}											
										}elsif (($alignCentralNas[$i-1] && $alignCentralNas[$i] eq $alignCentralNas[$i-1]) || ($alignCentralNas[$i+2] && $alignCentralNas[$i] eq $alignCentralNas[$i+2])) {
											$i++; # homopolymer stretch, set next position after the deletion stretch
										}else {
											$flag = 1;
											last;
										}
									}else {
										$flag = 1;
										last;
									}
								}elsif ($alignCentralNas[$i] eq '-') {	# single insertion
									if ($i == 0) {
										if ($alignRestNas[$i] eq $alignRestNas[$i+1]) {
											# homopolymer stretch, do nothing
										}else {
											$flag = 1;
											last;
										}										
									}elsif (($alignRestNas[$i-1] && $alignRestNas[$i] eq $alignRestNas[$i-1]) || ($alignRestNas[$i+1] && $alignRestNas[$i] eq $alignRestNas[$i+1])) {
										# homopolymer stretch, do nothing
									}else {
										$flag = 1;
										last;
									}
								}elsif ($alignRestNas[$i] eq '-') {	# single deletion
									if ($i == 0) {
										if ($alignCentralNas[$i] eq $alignCentralNas[$i+1]) {
											# homopolymer stretch, do nothing
										}else {
											$flag = 1;
											last;
										}										
									}elsif (($alignCentralNas[$i-1] && $alignCentralNas[$i] eq $alignCentralNas[$i-1]) || ($alignCentralNas[$i+1] && $alignCentralNas[$i] eq $alignCentralNas[$i+1])) {
										# homopolymer stretch, do nothing
									}else {
										$flag = 1;
										last;
									}
								}
							}
							
							unless ($flag) {
								$picked{$rest_name} = 1;
								push @cluster_names, $rest_name;
								$count++;
								$total_rest += $rest_duplicates;
								last;
							}						
						}
						if ($flag) {	# not picked, write distance into file
							print DIST "$central_seq\t$rest_seq\t$indels\t", join ("\t", @$etStrs), "\n";
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
		
		my $consensus = $nameSeq{$central_name};
		if ($total_rest > $central_duplicates) {
			my $localOutAlignFile = $localOutFile;
			$localOutAlignFile =~ s/\.fas/.afa/;
			system ("perl $scriptPath/alignRegion.pl -if $localOutFile -oa $localOutAlignFile -mt $match -mm $misMatch -gp $gapPenalty -uf");
			$consensus = utils::GetCons ($localOutAlignFile);
			$consensus =~ s/\-//g;
			if ($consensus ne $central_seq) {
				my $flag = 0;
				foreach my $name (@cluster_names) {
					my $seq = $nameSeq{$name};
					$seq =~ s/\-//g;
					if ($consensus eq $seq) {
						$flag = 1;
						last;
					}
				}
				if (!$flag) {
					$consensus = $nameSeq{$central_name};
				}
			}
		}
		foreach my $name (@cluster_names) {
			print OUTPUT ">$name\n$consensus\n";
		}
	}
}
close OUTPUT;
close DIST;
print "Done!\n";
