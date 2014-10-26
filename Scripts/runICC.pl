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
use Parallel::ForkManager;
use paths;
use utils;
use File::Basename;
use File::Copy;
use File::Path;
use Getopt::Long;
use Cwd;
use POSIX qw(strftime);
use Config;

my %option = (
	'od' => 'ICC_output',
	'cf' => 0.05,
	'cs' => 2,
	'proc' => 1,
	'mt' => 10,
	'mm' => -9,
	'gp' => -15,
	'u'  => 0.00013,
	'ar' => '',
	'h'  => '',
);

my $usage = "\nUsage: perl runICC.pl [-option value]

options:  
-od     output directory name (default: $option{'od'})
-cs     minimal cluster size as a cluster seed (defualt: $option{'cs'})
-cf     carry-forward cutoff (0 - 1, default: $option{'cf'}, carry forward whose frequency is fewer than cutoff will be corrected)
-u      overall mismatch rate per site (default: $option{'u'})
-proc	use the specified number of processors (default: 1)
-mt     match score (default: $option{mt})
-mm     mismatch score (default: $option{mm})
-gp     gap penalty (default: $option{gp})
-ar     assemble corrected reads in each window into full length reads
-h      usage help

";

GetOptions (\%option, 'od=s', 'cf=f', 'u=f', 'proc=i', 'mt=f', 'mm=f', 'gp=f', 'cs=i', 'ar', 'h');

my $outDir = $option{'od'} or die $usage;
my $cfCut = $option{'cf'};
my $minClusterSize = $option{'cs'};
my $u = $option{'u'};
my $proc = $option{'proc'};
my $match = $option{'mt'};
my $mismatch  = $option{'mm'};
my $gapPenalty = $option{'gp'};
my $assembleReads = $option{'ar'};
my $help = $option{'h'};
die $usage if $help;
my $inDir = getcwd();
my $scriptPath = $paths::scriptPath;

if ($proc < 1) {
	die "Number of processors must be greater or equal 1.\n";
}
my @workingDirs = ();
opendir DH, $inDir or die "couldn't open 1 $inDir\n";
while (my $subdir = readdir DH) {
	unless ($subdir =~ /^\./) {
		if ($subdir =~ /Region(\d+)\-(\d+)/) {
			push @workingDirs, $subdir;
		}
	}
}
closedir DH;
die "No region directorys exist. Working directory must be the parent of regions' directories.\n" if (!@workingDirs);

my $pref_name = '';
if ($inDir =~ /(.*)\/(.*?)$/) {
	$pref_name = $2;
}elsif ($inDir =~ /(.*)\\(.*?)$/) {
	$pref_name = $2;
}else {
	$pref_name = $inDir;
}

my $startTime = time();
my @sortDirs = sort by_number @workingDirs;

my $pm = Parallel::ForkManager->new($proc);
foreach my $subdir (@sortDirs) {
	my $pid = $pm->start and next;
	$subdir = $inDir.'/'.$subdir;	
	if (-d $subdir) {
		if ($subdir =~ /Region(\d+)\-(\d+)/i) {
			my $rgStart = $1;
			my $rgEnd = $2;
			my $indelCfCorrOutUniq = '';
			my $refSeq = '';
			my $refFile = my $logFile = my $nfFile = my $nhFile = my $nhfFile = my $anhfFile = my $ahFile = my $ahfFile = my $snvFile = '';
			my $fwdFlag = my $revFlag = 0;
			my $portionStart = my $portionEnd = 0;
			
			opendir SUB, $subdir or die "couldn't open $subdir\n";
			while (my $ssubdir = readdir SUB) {
				unless ($ssubdir =~ /^\./) {
					$ssubdir = $subdir.'/'.$ssubdir;
					if (-d $ssubdir && $ssubdir =~ /F_R_combo/) {						
						print "\n====== Processing $ssubdir ======\n";
						opendir SSUB, $ssubdir or die "couldn't open $ssubdir: $!\n";
						my $localOutDir = $ssubdir.'/'.$outDir;															
						if (-e $localOutDir) {
							rmtree($localOutDir);
						}
						mkdir $localOutDir;						
						chdir $localOutDir;
						while (my $file = readdir SSUB) {
							if ($file =~ /_combo\.fas$/) {
								my $inFile = $ssubdir.'/'.$file;
								my $uniqFile = $localOutDir.'/'.$file;
								$uniqFile =~ s/\.fas$/_U.fas/;								
								$refFile = $localOutDir.'/reference.fas';
								$logFile = $localOutDir.'/runICC.log';
								$nfFile = $localOutDir.'/NF.txt';
								$nhFile = $localOutDir.'/NH.txt';
								$nhfFile = $localOutDir.'/NHF.txt';
								$anhfFile = $localOutDir.'/ANHF.txt';
								$ahFile = $localOutDir.'/AH.txt';
								$ahfFile = $localOutDir.'/AHF.txt';
								$snvFile = $localOutDir.'/SNV.txt';
								my $refFlag = 0;
								open IN, $inFile or die "couldn't open $inFile: $!\n";
								open REF, ">$refFile" or die "couldn't open $refFile: $!\n";
								while (my $line = <IN>) {
									chomp $line;
									next if $line =~ /^\s*$/;
									if ($line =~ /^>Reference/i) {
										$refFlag = 1;
									}elsif ($refFlag) {
										$refSeq = $line;												
										print REF ">Reference\n$refSeq\n";
										last;
									}
								}
								close IN;
								close REF;
								# compress to unique reads
								system ("perl $scriptPath/uniqueReads.pl -if $inFile -of $uniqFile >$logFile");

								my $hmindelCorrOut = my $distOut = $uniqFile;
								$hmindelCorrOut =~ s/\.fas/_HIC.fas/;
								$distOut =~ s/\.fas/_dist.txt/;
								system ("perl $scriptPath/HIC.pl -if $uniqFile -of $hmindelCorrOut -od $distOut -cs $minClusterSize -mt $match -mm $mismatch -gp $gapPenalty >>$logFile");
							
								my $hmindelCorrOutUniq = $hmindelCorrOut;
								$hmindelCorrOutUniq =~ s/\.fas/_U.fas/;
								system ("perl $scriptPath/uniqueReads.pl -if $hmindelCorrOut -of $hmindelCorrOutUniq -uf >>$logFile");

								my $indelCorrOut = $hmindelCorrOutUniq;
								$indelCorrOut =~ s/_HIC_U\.fas/_IC.fas/;
								system ("perl $scriptPath/IC.pl -if $hmindelCorrOutUniq -id $distOut -of $indelCorrOut -cs $minClusterSize -mt $match -mm $mismatch -gp $gapPenalty >>$logFile");
								my $indelCorrOutUniq = $indelCorrOut;
								$indelCorrOutUniq =~ s/\.fas/_U.fas/;
								system ("perl $scriptPath/uniqueReads.pl -if $indelCorrOut -of $indelCorrOutUniq -uf >>$logFile");

								my $indelCfCorrOut = $indelCorrOutUniq;
								$indelCfCorrOut =~ s/_IC_U\.fas/_CC.fas/;
								$indelCfCorrOutUniq = $indelCfCorrOut;
								$indelCfCorrOutUniq =~ s/\.fas/_U.fas/;
								if ($cfCut == 0) {
									my $indelCfCorrOutUniqName = $indelCfCorrOutUniq.'.name';
									my $indelCorrOutUniqName = $indelCorrOutUniq.'.name';
									#copy($indelCorrOut, $indelCfCorrOut);
									copy($indelCorrOutUniq, $indelCfCorrOutUniq);
									open INNAME, $indelCorrOutUniqName or die "couldn't open $indelCorrOutUniqName: $!\n";
									open OUTNAME, ">", $indelCfCorrOutUniqName or die "couldn't open > $indelCfCorrOutUniqName: $!\n";
									while (my $line = <INNAME>) {
										chomp $line;
										next if $line =~ /^\s*$/;
										my ($name) = split /\s/, $line;
										print OUTNAME "$name $name\n";
									}
									close INNAME;
									close OUTNAME;
								}else {
									system ("perl $scriptPath/CC.pl -if $indelCorrOutUniq -of $indelCfCorrOut -cf $cfCut -mt 10 -mm -10 -gp -10 >>$logFile");								
									system ("perl $scriptPath/uniqueReads.pl -if $indelCfCorrOut -of $indelCfCorrOutUniq -uf >>$logFile");	
								}															
								unlink $distOut;			
							}
						}						
						closedir SSUB;
					}						
				}				
			}
			closedir SUB;
			
			my $indelCfCorrOutUniq_ref = $indelCfCorrOutUniq;
			$indelCfCorrOutUniq_ref =~ s/\.fas/_ref.fas/;
			open IN, $indelCfCorrOutUniq or die "couldn't open $indelCfCorrOutUniq: $!\n";
			open OUT, ">$indelCfCorrOutUniq_ref" or die "couldn't open $indelCfCorrOutUniq_ref: $!\n";
			print OUT ">Reference_1\n$refSeq\n";
			while (my $line = <IN>) {
				print OUT $line;
			}
			close IN;
			close OUT;
										
			# calculate Nt frequency
			my $afaFile = $indelCfCorrOutUniq_ref;
			$afaFile =~ s/\.fas/\.afa/;
			system ("perl $scriptPath/alignRegion.pl -if $indelCfCorrOutUniq_ref -oa $afaFile -mt $match -mm $mismatch -gp $gapPenalty -uf");		
			my $afaFreqFile = $afaFile;
			$afaFreqFile =~ s/\.afa/_freq\.txt/;
			system ("perl $scriptPath/ntFreq.pl -ia $afaFile -of $afaFreqFile -uf >>$logFile");
			
			# write frequencies in nucleotide and SNV frequency files with correct position
			open AFA, $afaFile or die "couldn't open $afaFile: $!\n";
			my $line = <AFA>;
			my $alignedRefSeq = <AFA>;
			my @alignedRefNas = split //, $alignedRefSeq;
			close AFA;
			open NF, ">$nfFile" or die "couldn't open $nfFile: $!\n";
			open NH, ">$nhFile" or die "could't open $nhFile: $!\n";
			open NHF, ">$nhfFile" or die "couldn't open $nhfFile: $!\n";
			open ANHF, ">$anhfFile" or die "couldn't open $anhfFile: $!\n";
			open AH, ">$ahFile" or die "couldn't open $ahFile: $!\n";
			open AHF, ">$ahfFile" or die "couldn't open $ahfFile: $!\n";
			open SNV, ">$snvFile" or die "couldn't open $snvFile: $!\n";			
			open FREQ, $afaFreqFile or die "couldn't open $afaFreqFile: $!\n";
			my $pos = $rgStart;
			my $afaFreqFlag = my $original_T_count = 0;
			while (my $line = <FREQ>) {
				chomp $line;
				next if ($line =~ /^\s*$/ || $line =~ /^Position/i);
				my @realPosFreqs = ();
				my ($idx, @rest) = split /\t/, $line;					
				++$afaFreqFlag;	
				if ($afaFreqFlag == 1) {
					$original_T_count = $rest[$#rest];
				}				
				if ($alignedRefNas[$idx-1] ne '-') {
					print NF "$pos\t";
					push @realPosFreqs, $pos, @rest;
					$pos++;						
				}else {
					print NF "\t";
					push @realPosFreqs, "", @rest;
				}
				print NF join ("\t", @rest), "\n";
				my $snvFreqs = utils::SNVcalling(\@realPosFreqs, $u);
				if (@$snvFreqs) {
					print SNV join("\t", @$snvFreqs), "\n";
				}
			}
			print NF "\n";
			print SNV "\n";
			close FREQ;
			
			# write to nucleotide haplotype fasta file
			my $poissonFlag = my $poissonCutReads = my $haploIdx = my $nt_dup = 0;
			my %ntSeqCount = my %aaSeqCount = ();
			my $poissonCut = utils::calculate_poisson($original_T_count, $u);
			open IN, $indelCfCorrOutUniq or die "couldn't open $indelCfCorrOutUniq: $!\n";			
			while (my $line = <IN>) {
				chomp $line;
				if ($line =~ /^>(.*?)_(\d+)$/) {
					$nt_dup = $2;
					if ($nt_dup >= $poissonCut) {
						$poissonFlag = 1;
						$haploIdx++;
						$poissonCutReads += $nt_dup;
						print NH ">Region".$rgStart."-".$rgEnd."_haplo".$haploIdx."_".$nt_dup."\n";
					}else {
						$poissonFlag = 0;
					}
				}elsif ($poissonFlag) {
					print NH "$line\n";
					$ntSeqCount{$line} = $nt_dup;
					my $idx = 0;
					my $aaSeq = '';
					while (length substr($line, $idx, 3) == 3) {
						my $codon = substr($line, $idx, 3);
						my $aa = '';
						if ($codon =~ /^[ACGT]{3}$/) {
							$aa = utils::translation($codon);
						}else {
							$aa = "X";							
						}
						$aaSeq .= $aa;
						$idx += 3;
					}
					$aaSeqCount{$aaSeq} += $nt_dup;
				}
			}
			close IN;
			
			# write to nucleotide haplotype frequency file
			foreach my $ntSeq (sort{$ntSeqCount{$b} <=> $ntSeqCount{$a}} keys %ntSeqCount) {
				my $freq = int ($ntSeqCount{$ntSeq} / $poissonCutReads * 1000000 + 0.5) / 1000000;
				print NHF "$ntSeq\t$rgStart\t$rgEnd\t$freq\t$ntSeqCount{$ntSeq}\n";
			}
			print NHF "\n";
			
			# write to aligned nucleotide haplotype frequency file
			my $refFlag = 0;
			open IN, $afaFile or die "couldn't open $afaFile: $!\n";
			while (my $line = <IN>) {
				chomp $line;
				if ($line =~ /^>Reference/) {
					$refFlag = 1;
				}elsif ($refFlag) {
					print ANHF "$line\t$rgStart\t$rgEnd\tRef\n";
					$refFlag = 0;
				}else {
					unless ($line =~ /^>/) {
						my $seq = $line;
						$seq =~ s/\-//g;
						if ($ntSeqCount{$seq}) {
							my $freq = int ($ntSeqCount{$seq} / $poissonCutReads * 1000000 + 0.5) / 1000000;
							print ANHF "$line\t$rgStart\t$rgEnd\t$freq\t$ntSeqCount{$seq}\n";
						}
					}
				}
			}
			print ANHF "+\n";
			close IN;
			
			# write to amino acid haplotype fasta and frequency file
			my $aaHaploIdx = 0;
			foreach my $aaSeq (sort{$aaSeqCount{$b} <=> $aaSeqCount{$a}} keys %aaSeqCount) {
				$aaHaploIdx++;
				my $freq = int ($aaSeqCount{$aaSeq} / $poissonCutReads * 1000000 + 0.5) / 1000000;
				print AHF "$aaSeq\t$rgStart\t$rgEnd\t$freq\t$aaSeqCount{$aaSeq}\n";
				print AH ">Region".$rgStart."-".$rgEnd."_haplo".$aaHaploIdx."_".$aaSeqCount{$aaSeq}."\n";
				print AH "$aaSeq\n";
			}
			print AHF "\n";
			close NF;
			close NH;
			close NHF;
			close ANHF;
			close AH;
			close AHF;
			close SNV;
		}			
	}
	$pm->finish; # Terminates the child process
}
$pm->wait_all_children;

chdir $inDir;

my $ntFreqFile = $pref_name."_nt_freq.txt";
my $ntHaploFile = $pref_name."_nt_haplotypes.fas";
my $ntHaploFreqFile = $pref_name."_nt_haplo_freq.txt";
my $alignedNtHaploFreqFile = $pref_name."_align_nt_haplo_freq.txt";
my $aaHaploFile = $pref_name."_aa_haplotypes.fas";
my $aaHaploFreqFile = $pref_name."_aa_haplo_freq.txt";
my $snvFreqFile = $pref_name."_SNV_freq.txt";
my $logFile = $pref_name."_ICC.log";

open NF, ">$ntFreqFile" or die "couldn't open $ntFreqFile: $!\n";
open NH, ">$ntHaploFile" or die "could't open $ntHaploFile: $!\n";
open NHF, ">$ntHaploFreqFile" or die "couldn't open $ntHaploFreqFile: $!\n";
open ANHF, ">$alignedNtHaploFreqFile" or die "couldn't open $alignedNtHaploFreqFile: $!\n";
open AH, ">$aaHaploFile" or die "couldn't open $aaHaploFile: $!\n";
open AHF, ">$aaHaploFreqFile" or die "couldn't open $aaHaploFreqFile: $!\n";
open SNV, ">$snvFreqFile" or die "couldn't open $snvFreqFile: $!\n";
open LOG, ">$logFile" or die "couldn't open $logFile: $!\n";
print NHF "Haplotype\tStart\tEnd\tFrequency\tReads\n";
print ANHF "Haplotype\tStart\tEnd\tFrequency\tReads\n";
print AHF "Haplotype\tStart\tEnd\tFrequency\tReads\n";
print NF "Position\tConsensus\t-\tA\tC\tG\tT\tCoverage\n";
print SNV "Position\tConsensus\t-\tA\tC\tG\tT\tCoverage\n";
print LOG "\n###########################################################################\n";
print LOG "# Input directory: $inDir\n";
print LOG "# Output directory: $outDir\n";
print LOG "# carryforward correction cut-off: $cfCut\n";
print LOG "# minimal cluster size as a cluster seed: $minClusterSize\n";
print LOG "# overall mismatch rate per site: $u\n";
print LOG "# number of processors: $proc\n";
print LOG "# match: $match; mismatch: $mismatch; gap penalty: $gapPenalty\n";
print LOG "###########################################################################\n\n";

foreach my $subdir (@sortDirs) {
	$subdir = $inDir.'/'.$subdir;
	print LOG "=== $subdir ===\n";		
	my $iccOutDir = $subdir.'/F_R_combo/'.$outDir;															
	if (-d $iccOutDir) {
		my $localLogFile = $iccOutDir.'/runICC.log';
		open LLOG, $localLogFile or die "couldn't open $localLogFile: $!\n";
		while (my $line = <LLOG>) {
			print LOG $line;
		}
		print LOG "\n";
		close LLOG;
		my $localNFFile = $iccOutDir.'/NF.txt';
		open LNF, $localNFFile or die "couldn't open $localNFFile: $!\n";
		while (my $line = <LNF>) {
			print NF $line;
		}
		close LNF;
		my $localNHFFile = $iccOutDir.'/NHF.txt';
		open LNHF, $localNHFFile or die "couldn't open $localNHFFile: $!\n";
		while (my $line = <LNHF>) {
			print NHF $line;
		}
		close LNHF;
		my $localANHFFile = $iccOutDir.'/ANHF.txt';
		open LANHF, $localANHFFile or die "couldn't open $localANHFFile: $!\n";
		while (my $line = <LANHF>) {
			print ANHF $line;
		}
		close LANHF;
		my $localNHFile = $iccOutDir.'/NH.txt';
		open LNH, $localNHFile or die "couldn't open $localNHFile: $!\n";
		while (my $line = <LNH>) {
			print NH $line;
		}
		close LNH;
		my $localAHFile = $iccOutDir.'/AH.txt';
		open LAH, $localAHFile or die "couldn't open $localAHFile: $!\n";
		while (my $line = <LAH>) {
			print AH $line;
		}
		close LAH;
		my $localAHFFile = $iccOutDir.'/AHF.txt';
		open LAHF, $localAHFFile or die "couldn't open $localAHFFile: $!\n";
		while (my $line = <LAHF>) {
			print AHF $line;
		}
		close LAHF;
		my $localSNVFile = $iccOutDir.'/SNV.txt';
		open LSNV, $localSNVFile or die "couldn't open $localSNVFile: $!\n";
		while (my $line = <LSNV>) {
			print SNV $line;
		}
		close LSNV;
	}else {
		print LOG "Couldn't find directory $iccOutDir.\n";
	}
}
close NF;
close NH;
close NHF;
close ANHF;
close AH;
close AHF;
close SNV;

if ($assembleReads) {
	my (@readnames, %readnameStatus, %readSeq, %corrreadSeq); 
	foreach my $subdir (@sortDirs) {
		if ($subdir =~ /Region(\d+)\-(\d+)/i) {			
			opendir SUB, $subdir or die "couldn't open $subdir, line 417\n";
			while (my $ssubdir = readdir SUB) {
				unless ($ssubdir =~ /^\./) {
					$ssubdir = $subdir.'/'.$ssubdir;
					if (-d $ssubdir && $ssubdir =~ /F_R_combo/) {						
						opendir SSUB, $ssubdir or die "couldn't open $ssubdir: $!\n";
						while (my $file = readdir SSUB) {
							if ($file =~ /_combo\.fas$/) {
								my $name = '';
								$file = "$ssubdir/$file";
								open COMBO, $file or die "couldn't open $file: $!\n";
								while (my $line = <COMBO>) {
									chomp $line;
									next if $line =~ /^\s*$/;
									if ($line =~ /^>(\S+)/) {
										$name = $1;
										if (!$readnameStatus{$name}) {
											$readnameStatus{$name} = 1;
											push @readnames, $name;
										}
									}else {
										$readSeq{$name} .= $line;
									}
								}
								close COMBO;
								last;
							}
						}
						closedir SSUB;
					}
				}
			}
			closedir SUB;
		}else {
			die "No region directory\n";
		}
		my @regionreadnames = ();				
		my $iccOutDir = $subdir.'/F_R_combo/'.$outDir;															
		if (-d $iccOutDir) {
			my $nameFile = my $hicnameFile = my $icnameFile = my $ccnameFile = my $ccreadFile = '';
			my %uniqReadnames = my %nameUniqname = my %uniqnameSeq = ();
			opendir ICC, $iccOutDir or die "couldn't open $iccOutDir: $!\n";
			while (my $file = readdir ICC) {
				if ($file =~ /combo_U.fas.name$/) {
					$nameFile = $file;
				}elsif ($file =~ /combo_U_HIC_U.fas.name$/) {
					$hicnameFile = $file;
				}elsif ($file =~ /combo_U_IC_U.fas.name$/) {
					$icnameFile = $file;
				}elsif ($file =~ /combo_U_CC_U.fas.name$/) {
					$ccnameFile = $file;
				}elsif ($file =~ /combo_U_CC_U.fas$/) {
					$ccreadFile = $file;
				}
			}
			closedir ICC;
			if ($nameFile) {
				$nameFile = $iccOutDir.'/'.$nameFile;
				open NAME, $nameFile or die "couldn't open $nameFile: $!\n";
				while (my $line = <NAME>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^(\S+)\s+(\S+)$/) {
						my $uniqName = $1;
						my $namestring = $2;
						my @names = split /,/, $namestring;						
						push @regionreadnames, @names;
						push @{$uniqReadnames{$uniqName}}, @names;
					}
				}
				close NAME;
			}else {
				die "Couldn't find combo_U.fas.name file\n";
			}
			if ($hicnameFile) {
				$hicnameFile = $iccOutDir.'/'.$hicnameFile;
				open NAME, $hicnameFile or die "couldn't open $hicnameFile: $!\n";
				while (my $line = <NAME>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^(\S+)\s+(\S+)$/) {
						my $uniqName = 'HIC_'.$1;
						my $namestring = $2;
						my @names = split /,/, $namestring;
						foreach my $name (@names) {
							push @{$uniqReadnames{$uniqName}}, @{$uniqReadnames{$name}};
							delete $uniqReadnames{$name};
						}
					}
				}
				close NAME;
			}else {
				die "Couldn't find combo_U_HIC_U.fas.name file\n";
			}
			if ($icnameFile) {
				$icnameFile = $iccOutDir.'/'.$icnameFile;
				open NAME, $icnameFile or die "couldn't open $icnameFile: $!\n";
				while (my $line = <NAME>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^(\S+)\s+(\S+)$/) {
						my $uniqName = 'IC_'.$1;
						my $namestring = $2;
						my @names = split /,/, $namestring;
						foreach my $name (@names) {
							$name = 'HIC_'.$name;
							push @{$uniqReadnames{$uniqName}}, @{$uniqReadnames{$name}};
							delete $uniqReadnames{$name};
						}
					}
				}
				close NAME;
			}else {
				die "Couldn't find combo_U_IC_U.fas.name file\n";
			}
			if ($ccnameFile) {
				$ccnameFile = $iccOutDir.'/'.$ccnameFile;
				my $outFile = $ccnameFile.'.out';
				open NAME, $ccnameFile or die "couldn't open $ccnameFile: $!\n";								
				open OUT, ">$outFile" or die "couldn't open $outFile: $!\n";
				while (my $line = <NAME>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^(\S+)\s+(\S+)$/) {
						my $origuniqname = $1;
						my $uniqName = 'CC_'.$1;
						my $namestring = $2;
						my @names = split /,/, $namestring;
						foreach my $name (@names) {
							$name = 'IC_'.$name;
							foreach my $readname (@{$uniqReadnames{$name}}) {
								$nameUniqname{$readname} = $origuniqname;
							}
							push @{$uniqReadnames{$uniqName}}, @{$uniqReadnames{$name}};
							delete $uniqReadnames{$name};
						}						
						print OUT $uniqName, " ", join(',', @{$uniqReadnames{$uniqName}}), "\n";
					}
				}
				close NAME;
			}else {
				die "Couldn't find combo_U_CC_U.fas.name file\n";
			}
			if ($ccreadFile) {
				my $name = '';
				$ccreadFile = $iccOutDir.'/'.$ccreadFile;
				open SEQ, $ccreadFile or die "couldn't open $ccreadFile: $!\n";
				while (my $line = <SEQ>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^>(\S+)/) {
						$name = $1;
					}else {
						$uniqnameSeq{$name} .= $line;
					}
				}
				close SEQ;
			}else {
				die "Couldn't find combo_U_CC_U.fas file\n";
			}
			foreach my $readname (@regionreadnames) {
				my $uniqname = $nameUniqname{$readname};
				my $seq = $uniqnameSeq{$uniqname};
				$corrreadSeq{$readname} .= $seq;				
			}
		}
	}
	my $aReadsb4Corr = $inDir.'/beforeCorr_assembled.fas';
	my $aReadsAfterCorr = $inDir.'/afterCorr_assembled.fas';
	open BC, ">$aReadsb4Corr" or die "couldn't open $aReadsb4Corr: $!\n";
	open AC, ">$aReadsAfterCorr" or die "couldn't open $aReadsAfterCorr: $!\n";
	foreach my $name (@readnames) {
		unless ($name =~ /Reference/i) {
			print BC ">$name\n$readSeq{$name}\n";
			print AC ">$name\n$corrreadSeq{$name}\n";
		}
		
	}
	close BC;
	close AC;
}

my $endTime = time();
my $duration = $endTime - $startTime;
my $hrs = int ($duration / 3600 * 100000 + 0.5) / 100000;
print LOG "\nAll done!\ntotal time duration: $hrs hour(s).\n\n";
close LOG;

sub by_number {
	$a =~ /Region(\d+)/;
	my $i = $1;
	$b =~ /Region(\d+)/;
	my $j = $1;
	$i <=> $j;
}
