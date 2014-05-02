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
use Cwd;
use lib "/home/wdeng/ICC/Scripts/lib";
use Parallel::ForkManager;

my %option = (
	'ref'    => '',
	'xmldir' => 'xmlOutput',
	'afc'    => 0.6,
	'rs'     => 0,
	're'     => 0,
	'proc'   => 1,
	'dlx'    => '',
	'h'      => '',
);

my $usage = "\nusage: perl retrieveRegion.pl [-option value]

options:  
-ref    input reference sequence used for blast in fasta format
-xmldir directory where BLASTN output XML files are (default: $option{'xmldir'})
-rs     region start position in reference sequence
-re     region end position in reference sequence
-afc    cutoff of read aligned fraction over the read length (default: $option{afc}, only above the cutoff the read is considered to be correctly aligned to reference)
-proc   number of processors (default: $option{proc})
-dlx    flag for deleting xml file after running the script (default: false)
-h      usage help
		
";

GetOptions (\%option, 'ref=s', 'xmldir=s', 'rs=i', 're=i', 'afc=f', 'proc=i', 'dlx', 'h');
my $refFile = $option{'ref'} or die $usage;
my $xmlDir = $option{'xmldir'};
my $startPos = $option{'rs'} or die $usage;
my $endPos = $option{'re'} or die $usage;
my $alignCut = $option{'afc'};
my $proc = $option{'proc'};
my $dlx = $option{'dlx'};
my $help = $option{'h'};
die $usage if $help;
my $inDir = getcwd();
my $pref_name = '';
if ($inDir =~ /(.*)\/(.*?)$/) {
	$pref_name = $2;
}elsif ($inDir =~ /(.*)\\(.*?)$/) {
	$pref_name = $2;
}else {
	$pref_name = $inDir;
}
print "xmlDir: $xmlDir\n";

my $refSeq = '';
my %readRgSeq = my %readRgSeqStatus = ();
my $queryCount = my $hitCount = my $passCutoffCount = my $notPassCutoffCount = my $rCutoffHitCount = my $notPassRcutHitCount = 0;
my $startTime = time();
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

my $pm = Parallel::ForkManager->new($proc);
# data structure retrieval and handling
$pm->run_on_finish ( 
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $retrievedData) = @_;
		# retrieve data structure from child
		if (defined($retrievedData)) {  # children are not forced to send anything
			#print "retrievedData: $retrievedData, ",scalar @$retrievedData,"\n";
			my $rgSeq = $retrievedData->[0];
			my $stats = $retrievedData->[1];
			foreach my $name (keys %$rgSeq) {
				if ($readRgSeqStatus{$name}) {
					die "duplicate key readRgSeqStatus{$name}\n";
				}else {
					$readRgSeqStatus{$name} = 1;
					$readRgSeq{$name} = $rgSeq->{$name};
				}
			}
			$queryCount += $stats->[0];
			$hitCount += $stats->[1];
			$passCutoffCount += $stats->[2];
			$notPassCutoffCount += $stats->[3];
			$rCutoffHitCount += $stats->[4];
			$notPassRcutHitCount += $stats->[5];
		}else {  # problems occuring during storage or retrieval will throw a warning
			print "No message received from child process $pid!\n";
		}
	}
);

opendir XMLDIR, $xmlDir or die "couldn't open $xmlDir: $!\n";
while (my $file = readdir XMLDIR) {		
	if ($file =~ /\.xml/) {
		my $pid = $pm->start and next;
		my $queryCount = my $queryLen = my $hitCount = my $hitFlag = my $hitStart = my $hitEnd = 0;
		my $qStart = my $qEnd = my $alignCutFlag = my $regionFlag = my $rCutoffHitCount = my $notPassRcutHitCount = my $frame = 0;
		my $passCutoffCount = my $notPassCutoffCount = 0;
		my $readName = my $alignedReadSeq = my $alignedRefSeq = '';
		my %regionSeq = ();
		my $xmlFile = "$xmlDir/$file";
		print "xmlFile: $xmlFile\n";
		open XML, $xmlFile or die "couldn't open $xmlFile: $!\n";
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
				}elsif ($frame == -1) {
					$readName = 'R_'.$readName;
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
		my @stats = ();
		my @data = ();
		push @stats, $queryCount, $hitCount, $passCutoffCount, $notPassCutoffCount, $rCutoffHitCount, $notPassRcutHitCount;
		push @data, \%regionSeq, \@stats;
		$pm->finish(0, \@data); # Terminates the child process
	}	
}
$pm->wait_all_children;
closedir XMLDIR;

my $dir = "Region_P".$startPos.'-'.$endPos;
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
my $fwdOutFile = $fwdDir.'/'.$pref_name.'_rg'.$startPos.'-'.$endPos.'_F.fas';
my $rvsOutFile = $rvsDir.'/'.$pref_name.'_rg'.$startPos.'-'.$endPos.'_R.fas';
my $comboOutFile = $comboDir.'/'.$pref_name.'_rg'.$startPos.'-'.$endPos.'_combo.fas';
open FWD, ">$fwdOutFile" or die "couldn't open $fwdOutFile: $!\n";
open RVS, ">$rvsOutFile" or die "couldn't open $rvsOutFile: $!\n";
open COMBO, ">$comboOutFile" or die "couldn't open $comboOutFile: $!\n";
print FWD ">Reference\n$refRgSeq\n";
print RVS ">Reference\n$refRgSeq\n";
print COMBO ">Reference\n$refRgSeq\n";
		
foreach my $name (keys %readRgSeq) {
	if ($name =~ /^F_/) {
		print FWD ">$name\n";
		print FWD $readRgSeq{$name},"\n";
	}elsif ($name =~ /^R_/) {
		print RVS ">$name\n";
		print RVS $readRgSeq{$name},"\n";
	}else {
		die "No direction info\n";
	}
	print COMBO ">$name\n";
	print COMBO $readRgSeq{$name},"\n";
}
close FWD;
close RVS;
close COMBO;

if ($dlx) {
	print "removing XML files ...\n";
	rmtree($xmlDir);
}
my $endTime = time();
my $duration = $endTime - $startTime;

print "Total queries: $queryCount, hits: $hitCount, hits pass align length percent cutoff of $alignCut: $passCutoffCount, not pass: $notPassCutoffCount\n";
print "hits cover the region of $startPos - $endPos: $rCutoffHitCount, not cover: $notPassRcutHitCount\n";
print "Time duration: $duration s.\n";


sub ReverseComplement {
	my $seq = shift;
	my $rcSeq = reverse $seq;
	$rcSeq =~ tr/ACGTacgt/TGCAtgca/;
	return $rcSeq;
}

