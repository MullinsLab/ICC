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
use File::Copy;
use File::Path;
use lib "/opt/home/wdeng/ICC/Scripts/lib";
use Parallel::ForkManager;

my %option = (
	'is'   => '',
	'iq'   => '',
	'os'   => '',
	'oq'   => '',
	'l'    => 100,
	'q'    => 25,
	'proc' => 1,
	'h'    => '', 
);

my $usage = "\nUsage: perl readQualFilter.pl [-option value]

options:  
-is     input reads fasta file
-iq     input quality file
-os     output reads fasta file
-oq     output quality file
-l      length cutoff (default: $option{'l'})
-q      average quality score cutoff (default: $option{'q'})
-proc   number of processors (default: $option{proc})
-h      show usage help

";

GetOptions (\%option, 'is=s', 'iq=s', 'os=s', 'oq=s', 'l=i','q=i', 'proc=i', 'h');

my $inFasta = $option{'is'} or die $usage;
my $inQual = $option{'iq'} or die $usage;
my $outFasta = $option{'os'} or die $usage;
my $outQual = $option{'oq'} or die $usage;
my $lenCut = $option{'l'};
my $qualCut = $option{'q'};
my $proc = $option{'proc'};
my $help = $option{'h'};
die $usage if $help;
unless ($lenCut =~ /^\d+$/) {
	die "length cutoff must be a positive integer number\n";
}
unless ($qualCut =~ /^\d+$/) {
	die "average quality cutoff must be a positive integer number\n";
}
unless ($proc =~ /^\d+$/ && $proc >= 1) {
	die "number of processors must be a positive integer number\n";
}
my $start = time();
my $readCount = my $totalSeq = my $totalQual = my $divideFastaFiles = my $divideQualFiles = 0;
my $belowLenCutCount = my $belowAllowNCount = my $belowQualCutCount = my $passFilterCount = my $failedCount = 0;
my $fastaInDir = 'fastaIn';
my $fastaOutDir = 'fastaOut';
my $qualInDir = 'qualIn';
my $qualOutDir = 'qualOut';
if (-e $fastaInDir) {
	rmtree($fastaInDir);
}
if (-e $fastaOutDir) {
	rmtree($fastaOutDir);
}
if (-e $qualInDir) {
	rmtree($qualInDir);
}
if (-e $qualOutDir) {
	rmtree($qualOutDir);
}
mkdir $fastaInDir;
mkdir $fastaOutDir;
mkdir $qualInDir;
mkdir $qualOutDir;

#if ($proc > 1) {	# split read file into small pieces
	# get total read number
	open IN, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <IN>) {
		if ($line =~ /^>/) {
			++$readCount;
		}
	}
	close IN;
	# divide reads into $proc small files
	my $files = $proc * 5;
	my $readsPerFile = int ($readCount / $files) + 1;
	$divideFastaFiles = DivideFile($readsPerFile, $inFasta, $fastaInDir, 'fas');
	$divideQualFiles = DivideFile($readsPerFile, $inQual, $qualInDir, 'qual');
	die "Something wrong with dividing files into smaller files\n" if ($divideFastaFiles != $divideQualFiles);
	
#}else {
#	copy($inFasta, "$fastaDir/");
#}

my $pm = Parallel::ForkManager->new($proc);
# data structure retrieval and handling
$pm->run_on_finish ( 
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $retrievedData) = @_;
		# retrieve data structure from child
		if (defined($retrievedData)) {  # children are not forced to send anything
			$totalSeq += $retrievedData->[0];
			$totalQual += $retrievedData->[1];
			$belowLenCutCount += $retrievedData->[2];
			$belowAllowNCount += $retrievedData->[3];
			$belowQualCutCount += $retrievedData->[4];
			$passFilterCount += $retrievedData->[5];
			$failedCount += $retrievedData->[6];
		}else {  # problems occuring during storage or retrieval will throw a warning
			print "No message received from child process $pid!\n";
		}
	}
);

for (my $i = 1; $i <= $divideFastaFiles; $i++) {
	my $pid = $pm->start and next;
	my $seq = my $seqName = "";
	my %failedName = ();
	my $seqCount = my $qualCount = 0;
	my $belowLenCutCount = my $belowAllowNCount = my $belowQualCutCount = my $passFilterCount = my $failedCount = 0;
	my $inFasta = "fastaIn/$i.fas";
	my $inQual = "qualIn/$i.qual";
	open INFASTA, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <INFASTA>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			if ($seq) {
				my $len = length $seq;
				if ($len < $lenCut) {
					$belowLenCutCount++;
					$failedName{$seqName} = 1;
				}elsif ($seq =~ /N/) {	# reads with N
					$belowAllowNCount++;
					$failedName{$seqName} = 1;					
				}
			}
			$seqName = $1;
			$seq = '';
			++$seqCount;
		}else {
			$seq .= uc $line;
		}
	}
	# last sequence
	if ($seq) {
		my $len = length $seq;
		if ($len < $lenCut) {
			$belowLenCutCount++;
			$failedName{$seqName} = 1;
		}elsif ($seq =~ /N/) {	# reads with N		
			$belowAllowNCount++;
			$failedName{$seqName} = 1;
		}
	}
	$seqName = $seq = '';
	close INFASTA;

	if ($qualCut > 0) {		
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
						$failedName{$seqName} = 1;
						$belowQualCutCount++;
					}
				}
				$seqName = $1;
				@quals = ();
				$totalQual++;
			}elsif (!$failedName{$seqName}) {
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
				$failedName{$seqName} = 1;
				$belowQualCutCount++;
			}
		}
		$seqName = '';
		close INQUAL;
	}

	my $outFasta = "fastaOut/$i.fas";
	open OUTFASTA, ">$outFasta" or die "couldn't open $outFasta: $!\n";
	open INFASTA, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <INFASTA>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			$seqName = $1;
			if (!$failedName{$seqName}) {
				print OUTFASTA ">$seqName\n";
				$passFilterCount++;
			}else {
				$failedCount++;
			}
		}elsif (!$failedName{$seqName}) {
			print OUTFASTA "$line\n";
		}
	}
	close OUTFASTA;
	close INFASTA;
	$seqName = '';

	my $outQual = "qualOut/$i.qual";
	open OUTQUAL, ">$outQual" or die "couldn't open $outQual: $!\n";
	open INQUAL, $inQual or die "couldn't open $inQual: $!\n";
	while (my $line = <INQUAL>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(\S+)/) {
			$seqName = $1;
			++$qualCount;
			if (!$failedName{$seqName}) {
				print OUTQUAL ">$seqName\n";
			}
		}elsif (!$failedName{$seqName}) {
			print OUTQUAL "$line\n";
		}
	}
	close OUTQUAL;
	close INQUAL;

	my @data = ();
	push @data, $seqCount, $qualCount, $belowLenCutCount, $belowAllowNCount, $belowQualCutCount, $passFilterCount, $failedCount;
	$pm->finish(0, \@data); # Terminates the child process	
}
$pm->wait_all_children;

open OUTFASTA, ">$outFasta" or die "couldn't open $outFasta: $!\n";
open OUTQUAL, ">$outQual" or die "couldn't open $outQual: $!\n";
for (my $i = 1; $i <= $divideFastaFiles; $i++) {
	my $inFasta = "fastaOut/$i.fas";
	my $inQual = "qualOut/$i.qual";
	open INFASTA, $inFasta or die "couldn't open $inFasta: $!\n";
	while (my $line = <INFASTA>) {
		print OUTFASTA $line;
	}
	close INFASTA;
	open INQUAL, $inQual or die "couldn't open $inQual: $!\n";
	while (my $line = <INQUAL>) {
		print OUTQUAL $line;
	}
	close INQUAL;
}
close OUTFASTA;
close OUTQUAL;

rmtree($fastaInDir);
rmtree($fastaOutDir);
rmtree($qualInDir);
rmtree($qualOutDir);

my $end = time();
my $duration = $end - $start;

print "total $totalSeq in seq file.\n$belowLenCutCount reads length shorter than $lenCut.\n";
print "$belowAllowNCount reads contain ambiguities.\n";
print "total $totalQual in qual file\n";
print "$belowQualCutCount reads average quality score below $qualCut.\n";
print "combined there are $passFilterCount meet the cutoff. $failedCount don't.\n";
print "Processor: $proc. Duration: $duration s. All done.\n";


sub DivideFile {
	my ($readsPerFile, $inFile, $dir, $ext) = @_;
	my $count = my $fileIdx = 0;
	my $fh;
	open IN, $inFile or die "couldn't open $inFile: $!\n";
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
				my $outFile = $dir.'/'.$fileIdx.'.'.$ext;
				open $fh, ">$outFile" or die "couldn't open $outFile: $!\n";
			}
		}
		print $fh $line,"\n";
	}
	close $fh;
	return $fileIdx;
}

