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
use lib "/opt/home/wdeng/ICC/Scripts/lib";
use Parallel::ForkManager;

my %option = (
	'ref'     => '',
	'rs'      => 1,
	're'      => 0,
	'afc'     => 0.6,
	'ws'      => 60,
	'ss'      => 60,
	'proc'    => 1,
	'dlx'     => '',
	'h'       => '',
);

my $usage = "\nusage: perl retrieveWindows.pl [-option value]

options:  
-ref    input reference sequence used for blast in fasta format
-rs     start position in reference sequence that the region start for retrieving windows (default: start at reference)
-re     end position in referece sequence that the region end for retrieving windows (default: end at reference)
-afc    cutoff of read aligned fraction over the read length (default: $option{afc}, only above the cutoff the read is considered to be correctly aligned to reference)
-ws     window size (default: $option{ws})
-ss     stride size (default: $option{ss})
-proc   number of processors (default: $option{proc})
-dlx    flag for deleting xml file after running the script (default: false)
-h      usage help
		
";

GetOptions (\%option, 'ref=s', 'rs=i', 're=i', 'afc=f', 'ws=i', 'ss=i', 'proc=i', 'dlx', 'h');
my $refFile = $option{'ref'} or die $usage;
my $wrStart = $option{'rs'};
my $wrEnd = $option{'re'};
my $alignCut = $option{'afc'};
my $wSize = $option{'ws'};
my $sSize = $option{'ss'};
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
my $xmlDir = 'xmlOutput';

my $refSeq = '';
my %readRgSeq = my %readRgSeqStatus = my %rgDirectionFlag = ();
my $queryCount = my $hitCount = my $passCutoffCount = my $notPassCutoffCount = my $wCutoffHitCount = my $notPassWcutHitCount = my $fullWindowReadCount = 0;
my $startTime = time();
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
			my $rgFlag = $retrievedData->[2];
			foreach my $i (keys %$rgSeq) {
				foreach my $name (keys %{$rgSeq->{$i}}) {
					if ($readRgSeqStatus{$i}{$name}) {
						die "duplicate key readRgSeqStatus{$i}{$name}\n";
					}else {
						$readRgSeqStatus{$i}{$name} = 1;
						$readRgSeq{$i}{$name} = $rgSeq->{$i}->{$name};
					}
				}
			}
			foreach my $i (keys %$rgFlag) {
				if ($rgFlag->{$i}->{F}) {
					$rgDirectionFlag{$i}{F} = $rgFlag->{$i}->{F};
				} 
				if ($rgFlag->{$i}->{R}) {
					$rgDirectionFlag{$i}{R} = $rgFlag->{$i}->{R};
				}
			}
			$queryCount += $stats->[0];
			$hitCount += $stats->[1];
			$passCutoffCount += $stats->[2];
			$notPassCutoffCount += $stats->[3];
			$wCutoffHitCount += $stats->[4];
			$notPassWcutHitCount += $stats->[5];
			$fullWindowReadCount += $stats->[6];
		}else {  # problems occuring during storage or retrieval will throw a warning
			print "No message received from child process $pid!\n";
		}
	}
);
opendir XMLDIR, $xmlDir or die "couldn't open $xmlDir: $!\n";
while (my $file = readdir XMLDIR) {		
	if ($file =~ /\.xml/) {
		my $pid = $pm->start and next;
		my $queryCount = my $queryLen = my $hitCount = my $hitFlag = my $start = my $end = 0;
		my $qStart = my $qEnd = my $alignCutFlag = my $windowFlag = my $wCutoffHitCount = my $notPassWcutHitCount = my $frame = 0;
		my $passCutoffCount = my $notPassCutoffCount = my $fullWindowReadCount = 0;
		my $readName = my $alignedReadSeq = my $alignedRefSeq = '';
		my (%readRegionSeq, %regionDirectionFlag);
		my $xmlFile = "$xmlDir/$file";
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
				}elsif ($frame == -1) {
					$readName = 'R_'.$readName;
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
				my $fullWindowReadFlag = 0;
				for (my $i = $wrStart; $i <= $wrEnd - $wSize + $sSize; $i += $sSize) {
					my $realWSize = $wSize;
					my $restLen = $wrEnd - $i + 1;
					if ($restLen < $wSize) {
						$realWSize = $restLen;
					}
					last if ($end < $i + $realWSize - 1);	# didn't cover last part
					next if ($start > $i);	# didn't cover biginning part
					if (!$fullWindowReadFlag) {
						$fullWindowReadFlag = 1;
						++$fullWindowReadCount;
					}
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
					$readRegionSeq{$i}{$readName} = substr($alignedReadSeq, $startIdx, $realLen);
					$readRegionSeq{$i}{$readName} =~ s/\-//g;
					
					if ($readName =~ /^F_/) {
						$regionDirectionFlag{$i}{'F'} = 1;
					}elsif ($readName =~ /^R_/) {
						$regionDirectionFlag{$i}{'R'} = 1;
					}
				}
			}
		}
		close XML;
		my @stats = ();
		my @data = ();
		push @stats, $queryCount, $hitCount, $passCutoffCount, $notPassCutoffCount, $wCutoffHitCount, $notPassWcutHitCount, $fullWindowReadCount;
		push @data, \%readRegionSeq, \@stats, \%regionDirectionFlag;
		$pm->finish(0, \@data); # Terminates the child process
	}	
}
$pm->wait_all_children;
closedir XMLDIR;

foreach my $regionStart (sort {$a <=> $b} keys %readRgSeq) {
	my $realWSize = $wSize;
	if ($regionStart + $wSize - 1 > $wrEnd) {
		$realWSize = $wrEnd - $regionStart + 1;
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
	my $fwdOutFile = $fwdDir.'/'.$pref_name.'_rg'.$regionStart.'-'.$regionEnd.'_F.fas';
	my $rvsOutFile = $rvsDir.'/'.$pref_name.'_rg'.$regionStart.'-'.$regionEnd.'_R.fas';
	my $comboOutFile = $comboDir.'/'.$pref_name.'_rg'.$regionStart.'-'.$regionEnd.'_combo.fas';
	if ($rgDirectionFlag{$regionStart}{F}) {
		open FWD, ">$fwdOutFile" or die "couldn't open $fwdOutFile: $!\n";
		print FWD ">Reference\n$refRgSeq\n";
	}
	if ($rgDirectionFlag{$regionStart}{R}) {
		open RVS, ">$rvsOutFile" or die "couldn't open $rvsOutFile: $!\n";
		print RVS ">Reference\n$refRgSeq\n";
	}		
	open COMBO, ">$comboOutFile" or die "couldn't open $comboOutFile: $!\n";	
	print COMBO ">Reference\n$refRgSeq\n";
			
	foreach my $name (keys %{$readRgSeq{$regionStart}}) {
		if ($name =~ /^F_/) {
			print FWD ">$name\n";
			print FWD $readRgSeq{$regionStart}{$name},"\n";
		}elsif ($name =~ /^R_/) {
			print RVS ">$name\n";
			print RVS $readRgSeq{$regionStart}{$name},"\n";
		}else {
			die "No direction info\n";
		}
		print COMBO ">$name\n";
		print COMBO $readRgSeq{$regionStart}{$name},"\n";
	}
	close FWD if ($rgDirectionFlag{$regionStart}{F});
	close RVS if ($rgDirectionFlag{$regionStart}{R});
	close COMBO;	
}
if ($dlx) {
	print "removing xml files ...\n";
	rmtree($xmlDir);
}
my $endTime = time();
my $duration = $endTime - $startTime;

print "Total queries: $queryCount, hits: $hitCount, hits pass align length percent cutoff of $alignCut: $passCutoffCount, not pass: $notPassCutoffCount\n";
print "hits pass window size of $wSize: $wCutoffHitCount, not pass: $notPassWcutHitCount, Read counts of covering at least one window: $fullWindowReadCount\n";
print "Time duration: $duration s.\n";



sub ReverseComplement {
	my $seq = shift;
	my $rcSeq = reverse $seq;
	$rcSeq =~ tr/ACGTacgt/TGCAtgca/;
	return $rcSeq;
}
