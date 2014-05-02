ICC 1.0
Contact: dengw@uw.edu

========
Overview
========

ICC is a software pipeline to correct sequencing errors such as indel and CAFIE errors
in Roche 454 pyrosequencing data, call single nucleotide variant and calculate 
haplotype frequencies.

===================
System Requirements
===================

ICC has been tested on systems running Linux, Mac OS X and MS Windows.

============
Installation
============

1. Download ICC_v1.0.zip to your preferred directory.
2. Unzip ICC_v1.0.zip.
3. Run "perl config.pl" in directory ICC_v1.0/Scripts.

Note: NCBI's BLAST+ packages (version 2.2.27+) for Linux, Mac and Windows have been 
pre-installed in ICC package. config.pl will automatically detect your operating system 
and configure the paths to run BLAST and perl scripts. If you move installed ICC_v1.0 
to other place, you need to run config.pl again before you execute ICC package.

=======================
How to run ICC package?
=======================

Starting with raw pyrosequencing reads fasta and quality files, ICC needs to run 
following four sequential steps:

1. Read quality filtering

Usage: perl /whereICCInstalled/Scripts/readQualFilter.pl [-option value]

Required parameters:
-is <file>            Input pyrosequencing reads fasta file
-iq <file>            Input pyrosequencing reads quality file
-os <file>            Output fasta file for reads passing quality filter
-oq <file>            Output quality file for reads passing quality filter

Other options:
-l <int>              Read length cutoff value. Read will be discard if it's length is
                      shorter than cutoff. [default: 100]
-q <int>              Read quality cutoff value. Read will be discard if it's average
                      quality is smaller than cutoff. [default: 25]
-h                    Usage help

2. Map filtered reads to reference by BLAST

Usage: perl /whereICCInstalled/Scripts/runBLAST.pl [-option value]

Required parameters:
-in <file>            Input quality filtered reads fasta file
-ref <file>           Input reference sequence fasta file
-out <file>           Output BLAST xml file

Other options:
-mt <int>             Match reward for BLAST algorithm. [default: 1]
-mm <int>             Mis-match penalty for BLAST algorithm. [default: -1]
-go <int>             Cost to open a gap for BLAST algorithm. [default: 1]  
-ge <int>             Cost to extend a gap for BLAST algorithm. [default: 2] 
-h                    Usage help         

3. Retrieve sequences in windows or individual regions

3.1. Retrieve consecutive windows' sequences

Usage: perl /whereICCInstalled/Scripts/retrieveWindows.pl [-option value]

Required parameters:
-xml <file>           Input BLAST xml file from step 2
-ref <file>           Input reference sequence fasta file

Other options:
-rs <int>             Start position in reference sequence to retrieve windows. 
                      [default: 1. i.e. beginning of the reference]
-re <int>             End position in reference sequence to retrieve windows. 
                      [default: 0. i.e. end of the reference]
-ws <int>             Window size. [default: 60]  
-ss <int>             Window stride size. [default: 60]  
-afc <int>            Alignment fraction Cutoff. Read will be considered to
                      specifically align to reference if the fraction of alingment is 
                      greater than cutoff.  [default: 0.6]  
-dlx                  Flag to delete xml file after running the script. [default: false] 
-h                    Usage help

Note: a serial sub-directories will be created in working directory with the name of
"Region_<num1>-<num2>". num1 is the start position of window in reference, num2 is 
the end position of window in reference. Each subdirectory contains retrieved read 
sequences covering the window. 

3.2. Retrieve sequences in a specific region

Usage: perl /whereICCInstalled/Scripts/retrieveRegion.pl [-option value]

Required parameters:
-xml <file>           Input BLAST xml file from step 2
-ref <file>           Input reference sequence fasta file
-rs <int>             Start position in reference sequence to retrieve region. 
-re <int>             End position in reference sequence to retrieve region. 

Other options:
-afc <float>          Alignment fraction Cutoff. Read will be considered to
                      specifically align to reference if the fraction of alignment is 
                      greater than cutoff.  [default: 0.6]  
-dlx                  Flag to delete xml file after running the script. [default: false]
-h                    Usage help

Note: a sub-directory will be created in working directory with the name of
"Region_<num1>-<num2>". num1 is the start position of the region in reference,
num2 is the end position of the region in reference. Subdirectory contains 
retrieved read sequences covering the region.  

4. Error correction, variant calling and profiling

Usage: perl /whereICCInstalled/Scripts/runICC.pl [-option value] >logFile

Options:
-od <file>            Output directory where ICC processing data and results will be 
                      stored. It will be created in the directory of 
                      yourWorkingDerectory/Region_<num1>-<num2>/F_R_combo/.
                      [default: ICC_output]
-cs <int>             Minimal cluster size as a cluster seed in three clustering steps
                      (homopolymer indel, indel and carry-forward errors only). 
                      Increasing the value will dramatically speed up the clustering.
                      [default: 2. i.e. clusters with only one read will not cluster
                      other reads]
-cf <float>           Carry-forward frequency cutoff. Carry-forward will be corrected
                      if it's frequency is lower than cutoff [default: 0.05]
-u <float>            Overall mismatch rate per site to approximate a Poisson 
                      distribution of error. [default: 0.00013]
-mt <int>             Match reward in multiple and pairwise alignment algorithm 
                      [default: 10]
-mm <int>             Mis-match penalty in multiple and pairwise alignment algorithm 
                      [default: -9]
-gp <int>             Gap penalty in multiple and pairwise alignment algorithm 
                      [default: -15]
-h                    Usage help

Note: You have to run runICC.pl in your working directory in which all subdirectories 
of "Region_<num1>-<num2>" are. It will output multiple files after the program finishes. 
logFile records the process of running the program. Other file names all begin with 
<NameOfYourWorkingDirectory>. _nt_freq.txt shows nucleotide frequencies at each position
across reference before variant calling using Poisson distribution model. _SNV_freq.txt
shows single nucleotide frequency by Poisson distribution. _nt_hyplotypes.fas is 
nucleotide hyplotype fasta file. _nt_hyplo_freq.txt is the file giving the frequency 
of each nucleotide hyplotype. _aa_hyplotypes.fas is amino acid hyplotype fasta file. 
_aa_hyplo_freq.txt is the file listing the frequency of each amino acid hyplotype. 

=======
Example
=======

The package includes example dataset you can test in the directory of ICC_v1.0/Example.

------------------------
How to run example data?
------------------------

1. Change working directory into ICC Example directory

Usage: cd /whereICCInstalled/Example

2. Example read quality filtering

Usage: perl /whereICCInstalled/Scripts/readQualFilter.pl -is exampleReads.fas -iq exampleReads.qual -os exampleReadsFilter.fas -oq exampleReadsFilter.qual

3. Map filtered reads to example reference by BLAST

Usage: perl /whereICCInstalled/Scripts/runBLAST.pl -in exampleReadsFilter.fas -ref exampleReference.fas -out exampleReadsFilterRef.xml

4. Retrieve consecutive windows' sequences

Usage: perl /whereICCInstalled/Scripts/retrieveWindows.pl -xml exampleReadsFilterRef.xml -ref exampleReference.fas

5. Error correction, variant calling and profiling

Usage: perl /whereICCInstalled/Scripts/runICC.pl  >exampleReadFilterRef.log

The following result files will be created: Example_aa_hyplo_freq.txt, Example_aa_hyplotypes.fas, 
Example_nt_freq.txt, Example_nt_hyplo_freq.txt, Example_nt_hyplotypes.fas, Example_SNV_freq.txt.
