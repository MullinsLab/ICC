package utils;

use strict;



sub SNVcalling {
	my $freqs = shift;
	my $u = shift;
	my %freq = ();
	my @snvFreqs = @$freqs;	
	my @nts = qw(- A C G T);	
	$freq{'-'} = $freqs->[2];
	$freq{'A'} = $freqs->[3];
	$freq{'C'} = $freqs->[4];
	$freq{'G'} = $freqs->[5];
	$freq{'T'} = $freqs->[6];
	my $readCount = $freqs->[7];

	my $restFreq = 0;
	foreach my $nt (@nts) {
		if ($nt ne $freqs->[1]) {
			$restFreq += $freq{$nt};
		}
	}
	if ($restFreq > 0) {
		my %ntCount = ();
		my $poissonCut = calculate_poisson($readCount, $u);
		my $freqCut = $poissonCut / $readCount;
		
		for (my $i = 0; $i < @nts; $i++) {			
			if ($freq{$nts[$i]} < $freqCut) {
				$freq{$nts[$i]} = 0;
			}
			$ntCount{$nts[$i]} = $readCount * $freq{$nts[$i]};				
		}
		my $finalCount = 0;
		foreach my $nt (@nts) {
			$finalCount += $ntCount{$nt}
		}
		$finalCount = int ($finalCount + 0.5);
		for (my $i = 0; $i < @nts; $i++) {
			$snvFreqs[$i+2] = int ($ntCount{$nts[$i]} / $finalCount * 1000000 + 0.5) / 1000000;
		}
		$snvFreqs[7] = $finalCount;
	}
	if ($snvFreqs[0] || $freq{'A'} > 0 || $freq{'C'} > 0 || $freq{'G'} > 0 || $freq{'T'} > 0) {
		
	}else {
		@snvFreqs = ();
	}
	return \@snvFreqs;
}

sub calculate_poisson {
	my $N = shift;
	my $u = shift;
	my $l = $N * $u;
	my $p = my $k = 0;
	my $P = my $d = 1;
	while ($P >= 0.001) {
		if ($k == 0 || $k == 1) {
			$d = 1;
		}else {
			$d = $k;
			for (my $i = $k-1; $i > 1; $i--) {
				$d = $d * $i;
			}
		}
		my $e = exp(-$l) * ($l**$k) / $d;
		$p += $e;
		$P = 1- $p;	
		++$k;
	}
	return $k;
}

sub translation {
	my $codon = shift;
	my %codons = (
		"ATT" => "I",
		"ATC" => "I",
		"ATA" => "I",
		"CTT" => "L",
		"CTC" => "L",
		"CTA" => "L",
		"CTG" => "L",
		"TTA" => "L",
		"TTG" => "L",
		"GTT" => "V",
		"GTC" => "V",
		"GTA" => "V",
		"GTG" => "V",
		"TTT" => "F",
		"TTC" => "F",
		"ATG" => "M",
		"TGT" => "C",
		"TGC" => "C",
		"GCT" => "A",
		"GCC" => "A",
		"GCA" => "A",
		"GCG" => "A",
		"GGT" => "G",
		"GGC" => "G",
		"GGA" => "G",
		"GGG" => "G",
		"CCT" => "P",
		"CCC" => "P",
		"CCA" => "P",
		"CCG" => "P",
		"ACT" => "T",
		"ACC" => "T",
		"ACA" => "T",
		"ACG" => "T",
		"TCT" => "S",
		"TCC" => "S",
		"TCA" => "S",
		"TCG" => "S",
		"AGT" => "S",
		"AGC" => "S",
		"TAT" => "Y",
		"TAC" => "Y",
		"TGG" => "W",
		"CAA" => "Q",
		"CAG" => "Q",
		"AAT" => "N",
		"AAC" => "N",
		"CAT" => "H",
		"CAC" => "H",
		"GAA" => "E",
		"GAG" => "E",
		"GAT" => "D",
		"GAC" => "D",
		"AAA" => "K",
		"AAG" => "K",
		"CGT" => "R",
		"CGC" => "R",
		"CGA" => "R",
		"CGG" => "R",
		"AGA" => "R",
		"AGG" => "R",
		"TAA" => '$',
		"TAG" => '$',
		"TGA" => '$',
	);	
	my $aa = $codons{$codon};
	return $aa;
}


sub GetCons {
	my $alignFile = shift;
	my @nas = ('-', 'A', 'G', 'C', 'T');
	my $duplicates = my $flag = my $alignLen = 0;
	my ($posNaCount, @consNas);
	open IN, $alignFile or die "couldn't open $alignFile: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^>(.*)_(\d+)$/) {
			$duplicates = $2;
		}else {
			if (!$flag) {
				$alignLen = length $line;
				$flag = 1;
			}else {
				my $len = length $line;
				if ($len != $alignLen) {
					die "not aligned\n";
				}
			}
			my @readNas = split //, uc $line;
			for (my $i = 0; $i < $alignLen; $i++) {
				$posNaCount->{$i}->{$readNas[$i]} += $duplicates;
			}
		}
	}
	for (my $i = 0; $i < $alignLen; $i++) {
		my $count = 0;
		my $consNa = '';
		foreach my $na (@nas) {
			if ($posNaCount->{$i}->{$na}) {
				if ($posNaCount->{$i}->{$na} > $count) {
					$consNa = $na;
					$count = $posNaCount->{$i}->{$na};
				}
			}
		}
		push @consNas, $consNa;
	}
	my $cons = join ('', @consNas);
	return $cons;
}



1;