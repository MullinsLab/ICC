package seqAlign;

use strict;


sub PairwiseAlign {
	my ($sub_central_seq, $sub_rest_seq, $sub_match, $sub_misMatch, $sub_gapPenalty) = @_;
	my @editTranscriptStr_array = ();
	my @centralNas = split //, $sub_central_seq;
	my @restNas = split //, $sub_rest_seq;
	my $centralLen = length $sub_central_seq;
	my $restLen = length $sub_rest_seq;
	my $value = CalculateSimilarity (\@centralNas, \@restNas, $centralLen, $restLen, $sub_match, $sub_misMatch, $sub_gapPenalty);		
	my $editTranscripts_preIns = Traceback_preIns (\@centralNas, \@restNas, $value, $centralLen, $restLen, $sub_match, $sub_misMatch, $sub_gapPenalty);
	my $editTranscriptStr_preIns = join ('', @$editTranscripts_preIns);
	my $MMCount = $editTranscriptStr_preIns =~ tr/R//;
	if (!$MMCount) {
		unless ($editTranscriptStr_preIns =~ /[ID]{3,}/) {	# three or more continuous indels, will not be corrected
			push @editTranscriptStr_array, $editTranscriptStr_preIns;
		}
		my $editTranscripts_preDel = Traceback_preDel (\@centralNas, \@restNas, $value, $centralLen, $restLen, $sub_match, $sub_misMatch, $sub_gapPenalty);
		my $editTranscriptStr_preDel = join ('', @$editTranscripts_preDel);
		my $MMCount2 = $editTranscriptStr_preDel =~ tr/R//;
		if (!$MMCount2) {
			unless ($editTranscriptStr_preIns eq $editTranscriptStr_preDel) {
				unless ($editTranscriptStr_preDel =~ /[ID]{3,}/) {	# three or more continuous indels, will not be corrected
					push @editTranscriptStr_array, $editTranscriptStr_preDel;
				}
			}
		}
	}
	return \@editTranscriptStr_array;
}

sub CalculateSimilarity {
	my ($consNasNoGapsRef, $readNasNoGapsRef, $consLen, $readLen, $match, $misMatch, $gapPenalty) = @_;
	my $value;
	for (my $i = 0; $i <= $consLen; $i++) {
		for (my $j = 0; $j <= $readLen; $j++) {
			if ($i == 0) {
				$value->[$i]->[$j] = $j * $gapPenalty;
			}elsif ($j == 0) {
				$value->[$i]->[$j] = $i * $gapPenalty;
			}else {
				my $upper = $value->[$i-1]->[$j] + $gapPenalty;
				my $left = $value->[$i]->[$j-1] + $gapPenalty;
				my $diagonal;
				if ($consNasNoGapsRef->[$i-1] eq $readNasNoGapsRef->[$j-1]) {
					$diagonal = $value->[$i-1]->[$j-1] + $match;
				}else {
					$diagonal = $value->[$i-1]->[$j-1] + $misMatch;
				}
				$value->[$i]->[$j] = Max ($upper, $left, $diagonal);
			}
		}
	}
	return $value;
}

sub Traceback_preIns {
	my ($consNasNoGaps, $readNasNoGaps, $value, $i, $j, $match, $misMatch, $gapPenalty) = @_;
	my @editTranscripts;	
	while ($i != 0 || $j != 0) {
		#print "i: $i, j: $j\n";
		my ($t, $editTranscript);
		if ($i != 0 && $j != 0) {
			if (@$consNasNoGaps[$i-1] eq @$readNasNoGaps[$j-1]) {
				$t = $match;
				$editTranscript = 'M';
			}else {
				$t = $misMatch;
				$editTranscript = 'R';
			}
			if ($value->[$i]->[$j] == $value->[$i-1]->[$j-1] + $t) {
				unshift @editTranscripts, $editTranscript;
				$i--;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i]->[$j-1] + $gapPenalty) {
				$editTranscript = 'I';
				unshift @editTranscripts, $editTranscript;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i-1]->[$j] + $gapPenalty) {
				$editTranscript = 'D';
				unshift @editTranscripts, $editTranscript;
				$i--;
			}
		}elsif ($j == 0) {
			$editTranscript = 'D';
			unshift @editTranscripts, $editTranscript;
			$i--;
		}elsif ($i == 0) {	
			$editTranscript = 'I';
			unshift @editTranscripts, $editTranscript;
			$j--;
		}
	}
	return \@editTranscripts;
}

sub Traceback_preDel {
	my ($consNasNoGaps, $readNasNoGaps, $value, $i, $j, $match, $misMatch, $gapPenalty) = @_;
	my @editTranscripts;	
	while ($i != 0 || $j != 0) {
		#print "i: $i, j: $j\n";
		my ($t, $editTranscript);
		if ($i != 0 && $j != 0) {
			if (@$consNasNoGaps[$i-1] eq @$readNasNoGaps[$j-1]) {
				$t = $match;
				$editTranscript = 'M';
			}else {
				$t = $misMatch;
				$editTranscript = 'R';
			}
			if ($value->[$i]->[$j] == $value->[$i-1]->[$j-1] + $t) {
				unshift @editTranscripts, $editTranscript;
				$i--;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i-1]->[$j] + $gapPenalty) {
				$editTranscript = 'D';
				unshift @editTranscripts, $editTranscript;
				$i--;
			}elsif ($value->[$i]->[$j] == $value->[$i]->[$j-1] + $gapPenalty) {
				$editTranscript = 'I';
				unshift @editTranscripts, $editTranscript;
				$j--;
			}
		}elsif ($j == 0) {
			$editTranscript = 'D';
			unshift @editTranscripts, $editTranscript;
			$i--;
		}elsif ($i == 0) {	
			$editTranscript = 'I';
			unshift @editTranscripts, $editTranscript;
			$j--;
		}
	}
	return \@editTranscripts;
}


sub Max {
	my @values = @_;
	my $flag = 0;
	my $max;
	foreach my $value (@values) {
		if (!$flag) {
			$max = $value;
			$flag++;
		}elsif ($value > $max) {
			$max = $value;
		}
	}
	return $max;
}


sub AlignRegion {
	my ($seqNamesRef, $nameSeqRef, $match, $misMatch, $gapPenalty, $gapMatch, $uniqFlag, $logFile) = @_;
	my @chars = qw(A C G T -);
	open LOG, ">$logFile" or die "couldn't open $logFile: $!\n";
	# reverse sequences so it will align homopolymers left to right after refinement
	foreach my $name (@$seqNamesRef) {
		$nameSeqRef->{$name} = reverse $nameSeqRef->{$name};
	}	
	my ($seqCount, $alignSeqs, $seqAlignseq, %nameAlignseq);	
	my ($uniqSeqs, $fullseqCount, $totalFullseqs) = GetUniqSeqs ($nameSeqRef, $uniqFlag);
	if (@$uniqSeqs) {
		print LOG "Total sequences: $totalFullseqs\n";
		if (@$uniqSeqs) {
			foreach my $uniqSeq (@$uniqSeqs) {
				print LOG "$uniqSeq: $fullseqCount->{$uniqSeq}\n";
			}
		}		
		if (@$uniqSeqs == 1) {	# only one unique full window sequence
			my $onlySeq = shift @$uniqSeqs;
			push @$alignSeqs, $onlySeq;
			$seqAlignseq->{$onlySeq} = $onlySeq;
			$seqCount->{$onlySeq} = $fullseqCount->{$onlySeq};
		}else {	# at least 2 unique full window sequences, first do pairwise alignment for the two most frequency sequences
			my $consSeq = shift @$uniqSeqs;
			my $readSeq = shift @$uniqSeqs;		
			my @consNasNoGaps = split //, $consSeq;
			my @readNasNoGaps = split //, $readSeq;
			my $consLen = length $consSeq;
			my $readLen = length $readSeq;
			my $value = CalculateSimilarity (\@consNasNoGaps, \@readNasNoGaps, $consLen, $readLen, $match, $misMatch, $gapPenalty);
			my $editTranscripts = Traceback_preIns (\@consNasNoGaps, \@readNasNoGaps, $value, $consLen, $readLen, $match, $misMatch, $gapPenalty);
			
			my $consAlignSeq = my $readAlignSeq = '';
			foreach my $et (@$editTranscripts) {
				if ($et eq 'I') {
					$consAlignSeq .= '-';
					$readAlignSeq .= shift @readNasNoGaps;
				}elsif ($et eq 'D') {
					$consAlignSeq .= shift @consNasNoGaps;
					$readAlignSeq .= '-';
				}else {
					$consAlignSeq .= shift @consNasNoGaps;
					$readAlignSeq .= shift @readNasNoGaps;
				}
			}
			$seqAlignseq->{$consSeq} = $consAlignSeq;
			$seqAlignseq->{$readSeq} = $readAlignSeq;
			$seqCount->{$consAlignSeq} = $fullseqCount->{$consSeq};
			$seqCount->{$readAlignSeq} = $fullseqCount->{$readSeq};
			push @$alignSeqs, $consAlignSeq, $readAlignSeq;
			die "problem in consNasNoGaps" if (@consNasNoGaps);
			die "problem in readNasNoGaps" if (@readNasNoGaps);						
			my $uniqcount = 2;
			# calculate initail profile
			my $profAlignLen = length $consAlignSeq;
			my $profSeqCount = $seqCount->{$consAlignSeq};
			my @consAlignNas = split //, $consAlignSeq;
			my $profileFreq;
			for (my $i = 1; $i <= $profAlignLen; $i++) {
				foreach my $char (@chars) {
					if ($consAlignNas[$i-1] eq $char) {
						$profileFreq->{$i}->{$char} = 1;
					}else {
						$profileFreq->{$i}->{$char} = 0;
					}
					
				}
			}			
			while (@$uniqSeqs) {	# at least 3 unique sequences, do profile alignment
				$uniqcount++;
				my $readSeq = shift @$uniqSeqs;
				my $readLen = length $readSeq;
				my @readNasNoGaps = split //, $readSeq;				
				my @profileSeqs = @$alignSeqs;				
				$profAlignLen = length $profileSeqs[0];
				my $newAlignSeq = $profileSeqs[$#profileSeqs];	# newly added aligned sequence in previous step, used to calculate the current profile
				RecalculateFreq($profileFreq, $newAlignSeq, $profSeqCount, \@chars, $seqCount, $profAlignLen);
				$profSeqCount += $seqCount->{$newAlignSeq};		
				foreach my $char (@chars) {
					print LOG "$char";
					for (my $i = 1; $i <= $profAlignLen; $i++) {
						my $freq = $profileFreq->{$i}->{$char} + 0.0000005;
						print LOG "\t";
						printf LOG ("%6f", $freq);
					}
					print LOG "\n";
				}
				print LOG "\n";		
				my $value = CalculateProfileSimilarity ($profileFreq, \@readNasNoGaps, \@chars, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch);
				my $editTranscripts = ProfileTraceback ($profileFreq, \@readNasNoGaps, \@chars, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPenalty, $gapMatch);				
				my $readAlignSeq = GetReadalignseq ($editTranscripts, \@readNasNoGaps);
				my $len = length $readAlignSeq;
				$seqCount->{$readAlignSeq} = $fullseqCount->{$readSeq};
				$seqAlignseq->{$readSeq} = $readAlignSeq;
								
				# re-align profile sequences and re-assign frequency if alignment length changes
				if (my $count = grep /D/, @$editTranscripts) {
					RealignProfile (\@profileSeqs, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq);
					my $pos = 0;
					my $realignProfFreq;
					for (my $i = 0; $i < @$editTranscripts; $i++) {
						$pos += 1;
						if ($editTranscripts->[$i] eq 'D') {
							foreach my $char (@chars) {
								if ($char eq '-') {
									$realignProfFreq->{$i+1}->{$char} = 1;
								}else {
									$realignProfFreq->{$i+1}->{$char} = 0;
								}								
							}
							$pos--;
						}else {
							foreach my $char (@chars) {
								$realignProfFreq->{$i+1}->{$char} = $profileFreq->{$pos}->{$char};							
							}
						}
					}
					$profileFreq = $realignProfFreq;
					undef $realignProfFreq;
				}
				push @$alignSeqs, $readAlignSeq;
			}			
		}
		print LOG "*** Aligned all sequences ***\n";
		print LOG join ("\n", @$alignSeqs), "\n";		
	}
	close LOG;	
	foreach my $seqName (@$seqNamesRef) {
		my $seq = $nameSeqRef->{$seqName};
		$nameAlignseq{$seqName} = $seqAlignseq->{$seq};
	}
	return \%nameAlignseq;
}

sub GetUniqSeqs {
	my $nameWindowseqRef = shift;
	my $uniqFlag = shift;
	my (%seqCount, @sortUniqSeqs);
	my $count = 0;
	foreach my $name (keys %$nameWindowseqRef) {
		my $alignSeq = my $pureSeq = $nameWindowseqRef->{$name};
		my $duplicates = 1;
		if ($uniqFlag) {
			$duplicates = GetDuplicates ($name);
		}
		$count += $duplicates;
		$pureSeq =~ s/\-//g;			
		if (!$seqCount{$pureSeq}) {
			$seqCount{$pureSeq} = 0;
		}
		$seqCount{$pureSeq} += $duplicates;
	}	
	foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		push @sortUniqSeqs, $seq;
	}
	return (\@sortUniqSeqs, \%seqCount, $count);
}


sub GetDuplicates {
	my $seqName = shift;
	my $duplicates;
	if ($seqName =~ /(\d+)$/) {
		$duplicates = $1;
	}else {
		die "there is no duplicate information at the end of sequence name\n";
	}
	return $duplicates;
}


sub RecalculateFreq {
	my ($sub_profileFreq, $sub_newAlignSeq, $sub_profSeqCount, $sub_chars, $sub_seqCount, $sub_profAlignLen) = @_;
	my @sub_newAlignNas = split //, $sub_newAlignSeq;
	my $totalSeqCount = $sub_profSeqCount + $sub_seqCount->{$sub_newAlignSeq};
	for (my $i = 1; $i <= $sub_profAlignLen; $i++) {
		foreach my $char (@$sub_chars) {
			my $charCount = $sub_profileFreq->{$i}->{$char} * $sub_profSeqCount;
			if ($sub_newAlignNas[$i-1] eq $char) {
				$charCount += $sub_seqCount->{$sub_newAlignSeq};
			}
			$sub_profileFreq->{$i}->{$char} = $charCount / $totalSeqCount;
		}
	}	
}


sub CalculateProfileSimilarity {
	my ($profileFreq, $readNasNoGapsRef, $charsRef, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch) = @_;
	my $value;
	for (my $i = 0; $i <= $readLen; $i++) {
		for (my $j = 0; $j <= $profAlignLen; $j++) {
			if ($i == 0 && $j == 0) {
				$value->[$i]->[$j] = 0;
			}elsif ($i == 0) {				
				my $current = 0;
				foreach my $char (@$charsRef) {
					if ($char eq '-') {
						$current += $gapMatch * $profileFreq->{$j}->{$char};
					}else {
						$current += $gapPenalty * $profileFreq->{$j}->{$char};
					}				
				}
				$value->[$i]->[$j] = $value->[$i]->[$j-1] + $current;				
			}elsif ($j == 0) {
				$value->[$i]->[$j] = $value->[$i-1]->[$j] + $gapPenalty;
			}else {
				my $gapCurrent = my $charCurrent = 0;
				foreach my $char (@$charsRef) {
					if ($char eq '-') {
						$charCurrent += $gapPenalty * $profileFreq->{$j}->{$char};
						$gapCurrent += $gapMatch * $profileFreq->{$j}->{$char};
					}else {
						$gapCurrent += $gapPenalty * $profileFreq->{$j}->{$char};
						if ($readNasNoGapsRef->[$i-1] eq $char) {
							$charCurrent += $match * $profileFreq->{$j}->{$char};
						}else {
							$charCurrent += $misMatch * $profileFreq->{$j}->{$char};
						}
					}
				}
				my $upper = $value->[$i-1]->[$j] + $gapPenalty;
				my $left = $value->[$i]->[$j-1] + $gapCurrent;
				my $diagonal = $value->[$i-1]->[$j-1] + $charCurrent;
				$value->[$i]->[$j] = Max ($upper, $left, $diagonal);
			}
		}
	}
	return $value;
}


sub ProfileTraceback {
	my ($profileFreq, $readNasNoGapsRef, $charsRef, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPanulty, $gapMatch) = @_;
	my @editTranscripts;
	my $idx = $profAlignLen;
	my $i = $readLen;
	my $j = $idx;
	while ($i != 0 || $j != 0) {
		my $editTranscript;
		my $gapCurrent = my $charCurrent = 0;
		if ($i != 0 && $j != 0) {
			foreach my $char (@$charsRef) {
				if ($char eq '-') {
					$charCurrent += $gapPanulty * $profileFreq->{$j}->{$char};
					$gapCurrent += $gapMatch * $profileFreq->{$j}->{$char};
				}else {
					$gapCurrent += $gapPanulty * $profileFreq->{$j}->{$char};
					if ($readNasNoGapsRef->[$i-1] eq $char) {
						$charCurrent += $match * $profileFreq->{$j}->{$char};
					}else {
						$charCurrent += $misMatch * $profileFreq->{$j}->{$char};
					}
				}
			}
			if ($value->[$i]->[$j] == $value->[$i-1]->[$j-1] + $charCurrent) {
				$editTranscript = 'M';
				unshift @editTranscripts, $editTranscript;
				$i--;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i]->[$j-1] + $gapCurrent) {
				$editTranscript = 'I';
				unshift @editTranscripts, $editTranscript;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i-1]->[$j] + $gapPanulty) {
				$editTranscript = 'D';
				unshift @editTranscripts, $editTranscript;
				$i--;
			}
		}elsif ($j == 0) {
			$editTranscript = 'D';
			unshift @editTranscripts, $editTranscript;
			$i--;
		}elsif ($i == 0) {	
			$editTranscript = 'I';
			unshift @editTranscripts, $editTranscript;
			$j--;
		}		
	}
	return \@editTranscripts;
}


sub GetReadalignseq {
	my ($editTranscripts, $readNasNoGapsRef) = @_;
	my $readAlignSeq = '';
	foreach my $et (@$editTranscripts) {
		if ($et eq 'I') {
			$readAlignSeq .= '-';
		}else {
			$readAlignSeq .= shift @$readNasNoGapsRef;
		}
	}
	return $readAlignSeq;
}


sub RealignProfile {
	my ($profileSeqsRef, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq) = @_;				
	for (my $i = 0; $i < @$profileSeqsRef; $i++) {
		my $profileSeq = $profileSeqsRef->[$i];
		my @profileSeqNas = split //, $profileSeq;
		my $alignSeq = '';
		foreach my $et (@$editTranscripts) {
			if ($et eq 'D') {
				$alignSeq .= '-';
			}else {
				$alignSeq .= shift @profileSeqNas;
			}
		}
		$alignSeqs->[$i] = $alignSeq;
		my $pureSeq = $profileSeq;
		$pureSeq =~ s/\-//g;		
		if ($seqCount->{$profileSeq}) {
			$seqCount->{$alignSeq} = $seqCount->{$profileSeq};
			$seqAlignseq->{$pureSeq} = $alignSeq;
		}	
	}
}


1;