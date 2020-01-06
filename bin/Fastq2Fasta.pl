#!/usr/bin/perl -w
# usage : Fastq2Fasta.pl fq
# note  : This converts a fastq file into a fasta file.

use strict;


# load fastq file and output in fasta format

open(IN,"<$ARGV[0]") || die "cannot open file $ARGV[0]!\n";

while(<IN>) {
    $_=">".substr($_,1);
    print $_;

    $_=<IN>;
    print $_;

    <IN>;
    <IN>;
}
close IN;
