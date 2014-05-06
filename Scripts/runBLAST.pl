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
use lib "/home/wdeng/ICC/Scripts/lib";
use paths;

my %option = (
	'in'  => '',
	'ref' => '',
	'out' => '',
	'mt'  => 1,
	'mm'  => -1,
	'go'  => 1,
	'ge'  => 2,
	'h'   => '',
);

my $usage = "\nUsage: perl runBLAST.pl [-option value]

options:  
-in      input reads fasta file
-out     BLAST output xml file
-ref     reference fasta file
-mt      match score (default: $option{'mt'})
-mm      mismatch score (default: $option{'mm'})
-go      cost to open a gap (default: $option{'go'})
-ge      cost to extend a gap (default: $option{'ge'})
-h       usage help

";

GetOptions (\%option, 'in=s', 'ref=s', 'out=s', 'bp=s', 'mt=i', 'mm=i','go=i', 'ge=i', 'h');

my $inFasta = $option{'in'} or die $usage;
my $refFasta = $option{'ref'} or die $usage;
my $outXML = $option{'out'} or die $usage;
my $reward = $option{'mt'};
my $penalty = $option{'mm'};
my $gapopen = $option{'go'};
my $gapextend = $option{'ge'};
my $help = $option{'h'};
die $usage if ($help);
my $blastPath = $paths::blastPath;

print "Making BLAST db ... ";
my $formatdbLog = $refFasta . ".log";
my $rv = 0;
$rv = system("$blastPath/makeblastdb -in $refFasta -dbtype nucl -logfile $formatdbLog");
unless ($rv == 0) {
	die "\nmakeblastdb failed: $rv\n";
}
print "done.\n";

print "BLASTing ... ";
$rv = system ("$blastPath/blastn -task blastn -db $refFasta -query $inFasta -out $outXML -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -dust no -outfmt 5");
unless ($rv == 0) {
	die "\nBLAST failed: $rv\n";
}
print "done.\n";
