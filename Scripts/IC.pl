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
use utils;
use paths;
use Getopt::Long;
use File::Path;

my %option = (
	'if' => '',
	'id' => '',
	'of' => '',
	'mt' => 10,
	'mm' => -9,
	'gp' => -15,
	'cs' => 2,
);

my $usage = "\nusage: perl IC.pl [-option value]

options:  
-if	input sequence fasta file
-id	input distance file
-of	output corrected sequence fasta file
-mt	match score (default: $option{'mt'})
-mm	mismatch score (default:  $option{'mm'})
-gp	gap penalty (default: $option{'gp'})
-cs	minimal cluster size as a cluster seed (defualt: $option{'cs'})
		
";

GetOptions (\%option, 'if=s', 'id=s', 'of=s', 'mt=f', 'mm=f', 'gp=f', 'cs=i');

my $inFile = $option{'if'} or die $usage;
my $inDist = $option{'id'} or die $usage;
my $outFile  = $option{'of'} or die $usage;
my $match = $option{'mt'};
my $misMatch = $option{'mm'};
my $gapPenalty = $option{'gp'};
my $minClusterSize = $option{'cs'};
my $scriptPath = $paths::scriptPath;
my $name = '';
my (@names, %nameSeq, %picked, $dist, $indel);
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

open DIST, $inDist or die "couldn't open $inDist: $!\n";
while (my $line = <DIST>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my ($first_seq, $snd_seq, $indels, @dists) = split /\t/, $line;
	$dist->{$first_seq}->{$snd_seq} = \@dists;
	$indel->{$first_seq}->{$snd_seq} = $indel->{$snd_seq}->{$first_seq} = $indels;
}
close DIST;

my $idx = my $count = 0;
my $outputDir = "./IC_otus";
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
					if ($dist->{$central_seq}->{$rest_seq} || $dist->{$rest_seq}->{$central_seq}) {
						
							$picked{$rest_name} = 1;
							push @cluster_names, $rest_name;
							my $rest_duplicates = 1;
							if ($rest_name =~ /_(\d+)$/) {
								$rest_duplicates = $1;
							}
							$total_rest += $rest_duplicates;
							$count++;
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
print "Done!\n";
