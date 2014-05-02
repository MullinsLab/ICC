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
use File::Basename;
use File::Path;
use Getopt::Long;

my %option = (
	'xml' => '',
	'ref' => '',
	'rs'  => 1,
	're'  => 0,
	'apc' => 0.6,
	'ws'  => 60,
	'ss'  => 60,
	'dlx' => '',
	'h'   => '',
);

my $usage = "\nusage: perl retrieveWindows.pl [-option value]

options:  
-xml    input blast output in xml format
-ref    input reference sequence used for blast in fasta format
-rs     start position in reference sequence that the region start for retrieving windows (default: start at reference)
-re     end position in referece sequence that the region end for retrieving windows (default: end at reference)
-apc    cutoff of read aligned part over the read length (default: $option{apc}, only above the cutoff the read is considered to be correctly aligned to reference)
-ws     window size (default: $option{ws})
-ss     stride size (default: $option{ss})
-dlx    flag for deleting xml file after running the script (default: false)
-h      usage help
		
";

GetOptions (\%option, 'xml=s', 'ref=s', 'rs=i', 're=i', 'apc=f', 'ws=i', 'ss=i', 'dlx', 'h');
my $xml = $option{'xml'} or die $usage;
my $refFile = $option{'ref'} or die $usage;
my $wrStart = $option{'rs'};
my $wrEnd = $option{'re'};
my $alignCut = $option{'apc'};
my $wSize = $option{'ws'};
my $sSize = $option{'ss'};
my $dlx = $option{'dlx'};
my $help = $option{'h'};
die $usage if $help;

my ($readName, $seq, %readAlignStart, %readAlignEnd, $referenceSeq, $readSeq, %alignLen, %refAlignNaLen, @alignReads);
my @refSeq;
my $fileName = basename($xml);
$fileName =~ s/\.(.*?)$//;
my $queryCount = my $queryLen = my $hitCount = my $fwdHitCount = my $rvsHitCount = my $hitFlag = my $flag = my $start = my $end = 0;
my $qStart = my $qEnd = my $alignCutFlag = my $windowFlag = my $wCutoffHitCount = my $notPassWcutHitCount = my $frame = 0;
my $passCutoffCount = my $notPassCutoffCount = my $lastWSize = 0;
my $refSeq = my $alignedReadSeq = my $alignedRefSeq = '';
my ($refRegionSeq, $refRegionStatus, $readRegionSeq, $readRegionStatus, $readRegionDup, %fwdregionFlag, %revregionFlag);
open REF, $refFile or die "couldn't open $refFile: $!\n";
while (my $line = <REF>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$line =~ s/\-//g;
		$refSeq .= $line;
	}
}
close REF;
my $refLen = length $refSeq;
print "refLen: $refLen\n";

if ($wrEnd == 0) {
	$wrEnd = $refLen;
}

open XML, $xml or die "couldn't open $xml: $!\n";
while (my $line = <XML>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /<Iteration_query-def>(.*)<\/Iteration_query-def>/) {
		$readName = $1;
		$queryCount++;
	}elsif ($line =~ /<Iteration_query-len>(\d+)<\/Iteration_query-len>/) {
		$queryLen = $1;
	}elsif ($line =~ /<Hit_num>1<\/Hit_num>/) {	# sometimes there are more than one hit, we just get the first one.
		$hitFlag = 1;
		$hitCount++;
	}elsif ($hitFlag && $line =~ /<Hsp_query-from>(\d+)<\/Hsp_query-from>/) {
		$qStart = $1;
	}elsif ($hitFlag && $line =~ /<Hsp_query-to>(\d+)<\/Hsp_query-to>/) {
		$qEnd = $1;
		my $qAlignLen = abs ($qEnd - $qStart) + 1;
		my $alignCoverage = $qAlignLen / $queryLen;
		if ($alignCoverage >= $alignCut) {
			$alignCutFlag = 1;
			$passCutoffCount++;
		}else {
			$hitFlag = $alignCutFlag = $windowFlag = 0;
			$notPassCutoffCount++;
			
		}
	}elsif ($alignCutFlag && $line =~ /<Hsp_hit-from>(\d+)<\/Hsp_hit-from>/) {
		$start = $1;
	}elsif ($alignCutFlag && $line =~ /<Hsp_hit-to>(\d+)<\/Hsp_hit-to>/) {
		$end = $1;
		my $alignLen = abs ($end - $start) + 1;
		if ($alignLen >= $wSize) {
			unless ($end > $start) {
				my $tmp = $start;
				$start = $end;
				$end = $tmp;
			}
			$windowFlag = 1;
			$wCutoffHitCount++;
		}else {
			$hitFlag = $alignCutFlag = $windowFlag = 0;
			$notPassWcutHitCount++;
		}
	}elsif ($windowFlag && $line =~ /<Hsp_hit-frame>(.*)<\/Hsp_hit-frame>/) {
		$frame = $1;		
		if ($frame == 1) {
			$readName = 'F_'.$readName;
			$fwdHitCount++;
		}elsif ($frame == -1) {
			$readName = 'R_'.$readName;
			$rvsHitCount++;
		}else {
			die "unrecognized frame of $frame\n";
		}				
	}elsif ($windowFlag && $line =~ /<Hsp_qseq>(.*)<\/Hsp_qseq>/) {
		$alignedReadSeq = $1;
	}elsif ($windowFlag && $line =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/) {	
		$alignedRefSeq = $1;
		if ($frame == -1) {	# reverse complement
			$alignedReadSeq = ReverseComplement ($alignedReadSeq);
			$alignedRefSeq = ReverseComplement ($alignedRefSeq);
		}
		my $alignedRefLen = length $alignedRefSeq;
		my @alignedRefNas = split //, $alignedRefSeq;
		$hitFlag = $alignCutFlag = $windowFlag = 0;

		for (my $i = $wrStart; $i <= $wrEnd - $wSize + $sSize; $i += $sSize) {
			my $realWSize = $wSize;
			my $restLen = $wrEnd - $i + 1;
			if ($restLen < $wSize) {
				$realWSize = $lastWSize = $restLen;
			}
			last if ($end < $i + $realWSize - 1);	# didn't cover last part
			next if ($start > $i);	# didn't cover biginning part
			my $preNaLen = my $startIdx = my $naLen = my $realLen = 0;
			my $preLen = $i - $start;
			while ($preNaLen < $preLen) {
				if ($alignedRefNas[$startIdx] =~ /[a-zA-Z]/) {
					$preNaLen++;
				}
				$startIdx++;
			}
			
			my $movingIdx = $startIdx;
			while ($naLen < $realWSize) {
				if ($alignedRefNas[$movingIdx] =~ /[a-zA-Z]/) {
					$naLen++;
				}
				$movingIdx++;
				$realLen++;
			}
			$readRegionSeq->{$i}->{$readName} = substr($alignedReadSeq, $startIdx, $realLen);
			$readRegionSeq->{$i}->{$readName} =~ s/\-//g;
			if ($readName =~ /^F_/) {
				$fwdregionFlag{$i} = 1;
			}elsif ($readName =~ /^R_/) {
				$revregionFlag{$i} = 1;
			}
		}
	}
}
close XML;
print "Total queries: $queryCount, hits: $hitCount, hits pass align length percent cutoff of $alignCut: $passCutoffCount, not pass: $notPassCutoffCount\n";
print "hits pass window size of $wSize: $wCutoffHitCount, not pass: $notPassWcutHitCount\n";

foreach my $regionStart (sort {$a <=> $b} keys %$readRegionSeq) {
	my $realWSize = $wSize;
	if ($regionStart + $wSize - 1 > $wrEnd) {
		$realWSize = $lastWSize;
	}
	my $refRgSeq = substr($refSeq, $regionStart - 1, $realWSize);
	my $regionEnd = $regionStart + $realWSize - 1;
	my $dir = "Region".$regionStart.'-'.$regionEnd;
	mkdir $dir unless (-e $dir);
	my $refFile = $dir.'/refSeq.fas';
	open REF, ">$refFile" or die "couldn't open $refFile: $!\n";
	print REF ">Reference\n$refRgSeq\n";
	
	my $fwdDir = $dir.'/F';	# directory storing forward reads
	my $rvsDir = $dir.'/R';	# directory storing reverse reads
	my $comboDir = $dir.'/F_R_combo';	# directory storing forward + reverse reads that contain forward/reverse information
	rmtree($fwdDir) if (-e $fwdDir);
	rmtree($rvsDir) if (-e $rvsDir);
	rmtree($comboDir) if (-e $comboDir);
	mkdir $fwdDir;
	mkdir $rvsDir;
	mkdir $comboDir;
	my $fwdOutFile = $fwdDir.'/'.$fileName.'_rg'.$regionStart.'-'.$regionEnd.'_F.fas';
	my $rvsOutFile = $rvsDir.'/'.$fileName.'_rg'.$regionStart.'-'.$regionEnd.'_R.fas';
	my $comboOutFile = $comboDir.'/'.$fileName.'_rg'.$regionStart.'-'.$regionEnd.'_combo.fas';
	if ($fwdregionFlag{$regionStart}) {
		open FWD, ">$fwdOutFile" or die "couldn't open $fwdOutFile: $!\n";
	}
	if ($revregionFlag{$regionStart}) {
		open RVS, ">$rvsOutFile" or die "couldn't open $rvsOutFile: $!\n";
	}		
	open COMBO, ">$comboOutFile" or die "couldn't open $comboOutFile: $!\n";	
	print FWD ">Reference\n$refRgSeq\n" if ($fwdregionFlag{$regionStart});
	print RVS ">Reference\n$refRgSeq\n" if ($revregionFlag{$regionStart});
	print COMBO ">Reference\n$refRgSeq\n";
			
	foreach my $name (keys %{$readRegionSeq->{$regionStart}}) {
		if ($name =~ /^F_/) {
			print FWD ">$name\n";
			print FWD $readRegionSeq->{$regionStart}->{$name},"\n";
		}elsif ($name =~ /^R_/) {
			print RVS ">$name\n";
			print RVS $readRegionSeq->{$regionStart}->{$name},"\n";
		}else {
			die "No direction info\n";
		}
		print COMBO ">$name\n";
		print COMBO $readRegionSeq->{$regionStart}->{$name},"\n";
	}
	close FWD if ($fwdregionFlag{$regionStart});
	close RVS if ($revregionFlag{$regionStart});
	close COMBO;	
}
if ($dlx) {
	print "removing $xml ...\n";
	unlink $xml;
}


sub ReverseComplement {
	my $seq = shift;
	my $rcSeq = reverse $seq;
	$rcSeq =~ tr/ACGTacgt/TGCAtgca/;
	return $rcSeq;
}