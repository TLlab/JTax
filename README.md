# JTax
Joining Illumina paired-end reads for classifying phylogenetic marker sequences

1. Introduction
-----------------------------------------------------------------------------------------------------
JTax is a tool designed to facilitate taxonomy classification for Illumina paired-end (PE) reads
of amplified marker genes (e.g., 16S rRNA) that cannot be merged into single reads either because
the high sequencing error rate at read ends prevents the merge or the paired reads are not long
enough to cover the whole amplicons. Given paired reads, JTax joins them in either a direct joining
or inside-out manner. The reference sequences are also rearranged correspondingly for classification.
JTax incorporates several classifiers, e.g., two word counting classifiers RDP classifier (rdp) and
SINTAX (sintax), and two top hit methods based on global alignment by USEARCH (top) and local
alignment by BLAST (btop).
-----------------------------------------------------------------------------------------------------


2. Installation
----------------------------------------------------------------------------------------------------
JTax is written in Perl and requires USEARCH, RDP Classifier, or BLAST depending on the selected
taxonomy classifiers. JTax has been tested on Linux and with USEARCH (v11), RDP classifier (v2.12),
and BLAST (v2.2.31+). Please see below for details of installation.
----------------------------------------------------------------------------------------------------

# JTax installation:
> wget https://github.com/TLlab/taxio/raw/master/jtax-0.1.tar.gz
> tar zxvf jtax-0.1.tar.gz
# add the folder jtax-0.1/bin to $PATH

# USEARCH installation:
# obtain USERACH (v11) from http://www.drive5.com/usearch/download.html
> chmod 755 usearch11.0.667_i86linux32
> ln -s usearch11.0.667_i86linux32 usearch
# add the folder containing usearch to $PATH

# RDP classifier installation:
> wget https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.12.zip
> unzip rdp_classifier_2.12.zip
# manually modify the path to RDP Classifier in jtax.pl, BuildReference.pl, and ClassifyJoinedRead.pl 

# BLAST installation:
> wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz
> tar zxvf ncbi-blast-2.2.31+-x64-linux.tar.gz
# add the folder ncbi-blast-2.2.31+/bin to $PATH


3. Running JTax
---------------------------------------------------------------------------------------------------------------------------------------------
usage  : jax.pl [options] ref.fa primer.fa read1.fq read2.fq
or     : jax.pl [options] ref.fa primer.fa read.fa/q

option : -method      <str> taxonomy classification method                                              [rdp*, nbc, sintax, ktop, top, btop]
       : -mrefid      <str> main reference ID                                                           []
       : -idmin       <flt> minimal identity of alignment between a reference and the main reference    [0.5; <=1]
       : -l1          <int> desired read 1 length                                                       [300; >=100]
       : -l2          <int> desired read 2 length                                                       [300; >=100]
       : -fr          <flt> reference length ratio to read length                                       [1.1; >1]
       : -trimp       <int> trim primer or not                                                          [1*,0]
       : -join        <str> joining method                                                              [io*,dj]
       : -minovlen    <int> minimal overlap length                                                      [16; >=10]
       : -pctid       <flt> minimal identity of overlap alignment                                       [0.75; <=1]
       : -nl          <int> length of padded Ns for direct joining                                      [8]
       : -out         <str> name of output folder                                                       []
       : -cpu         <int> number of threads                                                           [2]
       : -help              help
---------------------------------------------------------------------------------------------------------------------------------------------

JTax pipeline:
JTax takes a reference sequence file, a primer sequence file, and PE reads as input. Note that the input
reference shoud follow a format described below. Given the input, JTax first extracts amplicon segments
from the reference sequences and rearranges the amplicons in the direct joining or inside-out manner.
The rearranged references are then used for building a new database. JTax then joins the paired reads
accordingly either directly for taxonomy classification using the specified classifier.

To run JTax on the sample data, go to the test folder and issue the following commands, which creates
an output folder refio_test_rdp containing the inside-out references and inside-out reads, as well as
the classifications by RDP classifier. Check next section for descriptions of the output.

> jtax.pl ref.fa primer.fa test_1.fq test_2.fq

Example of classifying test inside-out reads against the inside-out references of RDP training data

> jtax.pl ../data/rdp16s.fa primer.fa test_1.fq test_2.fq

Example of applying the direct joining method and SINTAX for classification.

> jtax.pl -join dj -method sintax -cpu 8 ref.fa primer.fa test_1.fq test_2.fq

JTaxIO is designed to be modular, i.e., each step can be run separately as in the following example. 
This can be convenient when a user likes to classify additional data without rebulding the reference  
Results of the following commands should be identical to those generated above using RDP classifier
except for the slightly different confidence scores due to the random nature. 

> BuildReference.pl ref.fa primer.fa
> JoinPairedRead.pl test_1.fq test_2.fq > test.io.fq
> ClassifyJoinedRead.pl -out refio_testio.rdp refio test.io.fq


4. Input and output
-----------------------------------------------------------------------------------------------------
To run JTax successfully, the input data should be in a specific format as described below. 
-----------------------------------------------------------------------------------------------------

Input data

Reference sequence file:
JTax expects a reference file in FASTA format and the ID line should contain taxonomy information,
which should be formatted as ">accession;tax=d:domain,p:phylum,c:class,o:order,f:family,g:genus;". 
Please see the reference in the data or test folder for example.

Primer sequence file:
The primer file should also be in FASTA format and the ID of the forward and reverse primer should
contain the keyword FP and RP respectively. Please see the primer file in the test folder for
example. 

PE reads:
Input PE reads should be in FASTQ format (judged by the file extension fastq or fq). If single reads
are input, the reads can be in either FASTA or FASTQ format.

-----------------------------------------------------------------------------------------------------

Output data (using the sample data for examples)

rearranged reference: refio.fa or refdj.fa

joined reads: test.io.fq or test.dj.fq

classification results: refio_test.rdp or refdj_test.sintax

