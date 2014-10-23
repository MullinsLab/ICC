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
use File::Copy;
use File::Path;
use lib "/opt/home/wdeng/ICC/Scripts/lib";
use Parallel::ForkManager;
use paths;

my %option = (
	'in'   => '',
	'ref'  => '',
	'mt'   => 1,
	'mm'   => -1,
	'go'   => 1,
	'ge'   => 2,
	'proc' => 1,
	'h'    => '',
);

my $usage = "\nUsage: perl runBLAST.pl [-option value]

options:  
-in      input reads fasta file
-ref     reference fasta file
-mt      match score (default: $option{'mt'})
-mm      mismatch score (default: $option{'mm'})
-go      cost to open a gap (default: $option{'go'})
-ge      cost to extend a gap (default: $option{'ge'})
-proc    number of processors (default: $option{'proc'})
-h       usage help

";

GetOptions (\%option, 'in=s', 'ref=s', 'bp=s', 'mt=i', 'mm=i','go=i', 'ge=i', 'proc=i', 'h');

my $inFasta = $option{'in'} or die $usage;
my $refFasta = $option{'ref'} or die $usage;
my $reward = $option{'mt'};
my $penalty = $option{'mm'};
my $gapopen = $option{'go'};
my $gapextend = $option{'ge'};
my $proc = $option{'proc'};
my $help = $option{'h'};
die $usage if ($help);
my $blastPath = $paths::blastPath;
my $time1 = time();
print "Making BLAST db ... ";
my $formatdbLog = $refFasta . ".log";
my $rv = 0;
$rv = system("$blastPath/makeblastdb -in $refFasta -dbtype nucl -logfile $formatdbLog");
unless ($rv == 0) {
	die "\nmakeblastdb failed: $rv\n";
}
print "done.\n";

my $readCount = 0;
my $xmlDir = 'xmlOutput';
my $fastaDir = 'fastaOutput';
if (-e $xmlDir) {
	rmtree($xmlDir);
}
if (-e $fastaDir) {
	rmtree($fastaDir);
}
mkdir $xmlDir;
mkdir $fastaDir;
if ($proc > 1) {	# split read file into small pieces
	# get total read number
	open IN, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <IN>) {
		if ($line =~ /^>/) {
			++$readCount;
		}
	}
	close IN;
	# divide reads into $proc X 10 small files
	my $files = $proc * 10;
	my $readsPerFile = int ($readCount / $files) + 1;
	my $avgreads = $readCount / $files;
	my $count = my $fileIdx = 0;
	my $fh;
	open IN, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^>/) {
			++$count;
			if ($count > $readsPerFile) {
				$count = 1;
				close $fh;
			}
			if ($count == 1) {
				++$fileIdx;
				my $outFile = $fastaDir.'/'.$fileIdx.'.fas';
				open $fh, ">$outFile" or die "couldn't open $outFile: $!\n";
			}
		}
		print $fh $line,"\n";
	}
	close $fh;
}else {
	copy($inFasta, "$fastaDir/");
}
print "BLASTing ... ";
my $pm = Parallel::ForkManager->new($proc);
opendir FASTA, $fastaDir or die "couldn't open $fastaDir\n";
while (my $file = readdir FASTA) {	
	my $pid = $pm->start and next;
	unless ($file =~ /^\./) {
		my $inFile = $fastaDir.'/'.$file;
		my $outFile = $xmlDir.'/'.$file.'.xml';
		#print "$blastPath/blastn -task blastn -db $refFasta -query $inFile -out $outFile -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -dust no -outfmt 5\n";
		$rv = system ("$blastPath/blastn -task blastn -db $refFasta -query $inFile -out $outFile -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -dust no -outfmt 5");		
		unless ($rv == 0) {
			die "\nBLAST failed: $rv\n";
		}
	}
	$pm->finish; # Terminates the child process
}
$pm->wait_all_children;
closedir FASTA;
if (-e $fastaDir) {
	rmtree($fastaDir);
}
my $time2 = time();
my $duration = $time2 - $time1;
print "\nproc: $proc\n";
print "duration: $duration s. All done!\n";
