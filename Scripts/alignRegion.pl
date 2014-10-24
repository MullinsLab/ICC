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
use lib "/home/wdeng/ICC/Scripts/lib";
use seqAlign;
use Getopt::Long;

my %option = (
	'if' => '',
	'oa' => '',
	'mt' => 10,
	'mm' => -9,
	'gp' => -15,
	'uf' => ''
);

my $usage = "\nusage: perl alignRegion.pl [-option value]

options:  
-if		input alignment sequence fasta file
-oa		output refined alignment file
-mt		match score (default: $option{mt})
-mm		mismatch score (default: $option{mm})
-gp		gap penalty (default: $option{gp})
-uf		unique sequence flag (default: false). indicating if the sequences in alignment are unique,
		if true, the duplicate info must be at the end of sequence name
		
";

GetOptions (\%option, 'if=s', 'oa=s', 'mt=f', 'mm=f', 'gp=f', 'uf');

my $inFile = $option{'if'} or die $usage;
my $alignFile = $option{'oa'} or die $usage;
my $match = $option{'mt'};
my $misMatch = $option{'mm'};
my $gapPenalty = $option{'gp'};
my $uniqFlag = $option{'uf'};
my $gapMatch = 0;

my $logFile = $alignFile;
$logFile =~ s/\.\w+$/\.log/;
my (@seqNames, %nameSeq);
my $name = '';
open IN, $inFile or die "couldn't open $inFile: $!\n"; 
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		push @seqNames, $name;
	}else {
		$line =~ s/\-//g;
		$nameSeq{$name} .= uc $line;
	}
}
close IN;
	
# align sequences in the region
my $nameAlignseq = seqAlign::AlignRegion (\@seqNames, \%nameSeq, $match, $misMatch, $gapPenalty, $gapMatch, $uniqFlag, $logFile);
open ALGN, ">$alignFile" or die "couldn't open $alignFile: $!\n";
foreach my $seqName (@seqNames) {
	print ALGN ">$seqName\n";
	my $reverseBackRead = reverse $nameAlignseq->{$seqName};
	print ALGN $reverseBackRead,"\n";
}
close ALGN;
