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
use Cwd;
use Config;
use File::Copy;

my $scriptPath = getcwd();
my $os = $Config{osname};
$scriptPath =~ s/\/$//;
my $blastPath = $scriptPath;
$blastPath =~ s/Scripts/Blast/;
if ($os =~ /linux/i || $os =~ /unix/i) {
	$os = "Linux";
}elsif ($os =~ /darwin/i) {
	$os = "Mac";
}elsif ($os =~ /win/i) {
	$os = "Windows";
}else {
	die "Couldn't detect your operating system. Please contact dengw\"@\"uw.edu.\n";
}
$blastPath .= "/$os/bin";
my $lib = $scriptPath."/lib";
open PM, ">lib/paths.pm" or die "couldn't open paths.pm: $!\n";
print PM "package paths;\n\nuse strict;\n\n";
print PM "our \$scriptPath = \"$scriptPath\";\n";
print PM "our \$blastPath = \"$blastPath\";\n\n";
print PM "1;\n";
close PM;

my @pls = qw (runBLAST.pl runICC.pl alignRegion.pl HIC.pl IC.pl CC.pl);
foreach my $pl (@pls) {
	my $tmp = $pl;
	$tmp =~ s/\.pl$/_tmp.pl/;
	open IN, $pl or die "couldn't open $pl: $!\n";
	open OUT, ">$tmp" or die "couldn't open $tmp: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^use lib/) {
			print OUT "use lib \"$scriptPath/lib\";\n";
		}else {
			print OUT "$line\n";
		}
	}
	close IN;
	close OUT;
	move ($tmp, $pl);
	chmod 0755, $pl;
}
#print "os: $os\n";
#print "scriptPath: $scriptPath\n";
#print "blastPath: $blastPath\n";

