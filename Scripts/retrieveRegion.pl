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
	'afc' => 0.6,
	'rs'  => 0,
	're'  => 0,
	'dlx' => '',
	'h'   => '',
);

my $usage = "\nusage: perl retrieveRegion.pl [-option value]

options:  
-xml    input blast output in xml format
-ref    input reference sequence used for blast in fasta format
-rs     region start position in reference sequence
-re     region end position in reference sequence
-afc    cutoff of read aligned fraction over the read length (default: $option{afc}, only above the cutoff the read is considered to be correctly aligned to reference)
-dlx    flag for deleting xml file after running the script (default: false)
-h      usage help
		
";

GetOptions (\%option, 'xml=s', 'ref=s', 'rs=i', 're=i', 'afc=f', 'dlx', 'h');
my $xml = $option{'xml'} or die $usage;
my $refFile = $option{'ref'} or die $usage;
my $startPos = $option{'rs'} or die $usage;
my $endPos = $option{'re'} or die $usage;
my $alignCut = $option{'afc'};
my $dlx = $option{'dlx'};
my $help = $option{'h'};
die $usage if $help;

my ($readName, $seq, %readAlignStart, %readAlignEnd, $referenceSeq, $readSeq, %alignLen, %refAlignNaLen, @alignReads);
my @refSeq;
my $fileName = basename($xml);
$fileName =~ s/\.(.*?)$//;
my $queryCount = my $queryLen = my $hitCount = my $fwdHitCount = my $rvsHitCount = my $hitFlag = my $flag = my $hitStart = my $hitEnd = 0;
my $qStart = my $qEnd = my $alignCutFlag = my $regionFlag = my $rCutoffHitCount = my $notPassRcutHitCount = my $frame = 0;
my $passCutoffCount = my $notPassCutoffCount = 0;
my $refSeq = my $alignedReadSeq = my $alignedRefSeq = '';
my ($refRegionSeq, $refRegionStatus, $readRegionSeq, $readRegionStatus, $readRegionDup);
my %regionSeq;
open REF, $refFile or die "couldn't open $refFile: $!\n";
while (my $line = <REF>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$refSeq .= $line;
	}
}
close REF;
$refSeq =~ s/\-//g;
my $refRgSeq = substr($refSeq, $startPos-1, $endPos-$startPos+1);
my $refLen = length $refSeq;
print "refLen: $refLen\n";

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
			$hitFlag = $alignCutFlag = $regionFlag = 0;
			$notPassCutoffCount++;
			
		}
	}elsif ($hitFlag && $line =~ /<Hsp_hit-from>(\d+)<\/Hsp_hit-from>/) {
		$hitStart = $1;
	}elsif ($hitFlag && $line =~ /<Hsp_hit-to>(\d+)<\/Hsp_hit-to>/) {
		$hitEnd = $1;
		unless ($hitEnd > $hitStart) {
			my $tmp = $hitStart;
			$hitStart = $hitEnd;
			$hitEnd = $tmp;
		}
		if ($hitStart <= $startPos && $hitEnd >= $endPos) {
			$regionFlag = 1;
			$rCutoffHitCount++;
		}else {
			$hitFlag = $alignCutFlag = $regionFlag = 0;
			$notPassRcutHitCount++;
		}
	}elsif ($regionFlag && $line =~ /<Hsp_hit-frame>(.*)<\/Hsp_hit-frame>/) {
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
	}elsif ($regionFlag && $line =~ /<Hsp_qseq>(.*)<\/Hsp_qseq>/) {
		$alignedReadSeq = $1;
		
	}elsif ($regionFlag && $line =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/) {	
		$alignedRefSeq = $1;
		if ($frame == -1) {	# reverse complement
			$alignedReadSeq = ReverseComplement ($alignedReadSeq);
			$alignedRefSeq = ReverseComplement ($alignedRefSeq);
		}
		
		my $alignedRefLen = length $alignedRefSeq;
		my @alignedRefNas = split //, $alignedRefSeq;
		$hitFlag = $alignCutFlag = $regionFlag = 0;
		
		my $idx = my $startIdx = my $stopIdx = 0;
		for (my $i = 0; $i < $alignedRefLen; $i++) {
			if ($alignedRefNas[$i] =~ /[A-Za-z]/) {
				$idx++;
				if ($idx == $startPos-$hitStart+1) {
					$startIdx = $i;
				}elsif ($idx == $endPos-$hitStart+1) {
					$stopIdx = $i;
					last;
				}
			}
		}
		$regionSeq{$readName} = substr($alignedReadSeq, $startIdx, $stopIdx-$startIdx+1);
		$regionSeq{$readName} =~ s/\-//g;
	}
}
close XML;

print "Total queries: $queryCount, hits: $hitCount, hits pass align length percent cutoff of $alignCut: $passCutoffCount, not pass: $notPassCutoffCount\n";
print "hits cover the region of $startPos - $endPos: $rCutoffHitCount, not cover: $notPassRcutHitCount\n";

my $dir = "Region".$startPos.'-'.$endPos;
mkdir $dir unless (-e $dir);
my $regionRefFile = $dir.'/refSeq.fas';
open REF, ">$regionRefFile" or die "couldn't open $regionRefFile: $!\n";
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
my $fwdOutFile = $fwdDir.'/'.$fileName.'_rg'.$startPos.'-'.$endPos.'_F.fas';
my $rvsOutFile = $rvsDir.'/'.$fileName.'_rg'.$startPos.'-'.$endPos.'_R.fas';
my $comboOutFile = $comboDir.'/'.$fileName.'_rg'.$startPos.'-'.$endPos.'_combo.fas';
open FWD, ">$fwdOutFile" or die "couldn't open $fwdOutFile: $!\n";
open RVS, ">$rvsOutFile" or die "couldn't open $rvsOutFile: $!\n";
open COMBO, ">$comboOutFile" or die "couldn't open $comboOutFile: $!\n";
print FWD ">Reference\n$refRgSeq\n";
print RVS ">Reference\n$refRgSeq\n";
print COMBO ">Reference\n$refRgSeq\n";
		
foreach my $name (keys %regionSeq) {
	if ($name =~ /^F_/) {
		print FWD ">$name\n";
		print FWD $regionSeq{$name},"\n";
	}elsif ($name =~ /^R_/) {
		print RVS ">$name\n";
		print RVS $regionSeq{$name},"\n";
	}else {
		die "No direction info\n";
	}
	print COMBO ">$name\n";
	print COMBO $regionSeq{$name},"\n";
}
close FWD;
close RVS;
close COMBO;

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

