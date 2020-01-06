#!/usr/bin/perl -w
# usage : SplitFastq.pl fa/q n
# note  : This splits a fasta/q file into n files with about equal number of sequences.

use strict;
use POSIX qw(ceil);


# judge input file

my $f = $ARGV[0]=~/a$/ ? "fa" : "fq";


# link file if no splitting is required

if($ARGV[1]==1) {
    `ln -s $ARGV[0] $ARGV[0].1.$f`;
    exit 0;
}


# count number of sequences

my $nos;

if( $f eq "fa" ) {
    my $grep = `grep -c \\> $ARGV[0]`;
    ($nos) = $grep =~ /(\d+)/;

} else {
    my $wc = `wc -l $ARGV[0]`;
    ($nos) = $wc =~ /^(\d+)/; $nos/=4;
}
$nos = ceil( $nos / $ARGV[1] );


# load sequences and output

my $k = 1;
open OUT, ">$ARGV[0].$k.$f" || die "open $ARGV[0].$k.$f: $!\n";

my $i = 0;
open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

if( $f eq "fa" ) {
    while(<IN>) {
	$i++ if /^>/;
	if( $i > $nos ) {
	    close OUT;
	    $k++;
	    open OUT, ">$ARGV[0].$k.fa" || die "open $ARGV[0].$k.fa: $!\n";
	    $i = 1;
	}
	print OUT $_;
    }
} else {
    while(<IN>) {
	$i++;
	if( $i > $nos ) {
	    close OUT;
	    $k++;
	    open OUT, ">$ARGV[0].$k.$f" || die "open $ARGV[0].$k.$f: $!\n";
	    $i = 1;
	}
	print OUT $_;
	$_ = <IN>; print OUT $_;
	$_ = <IN>; print OUT $_;
	$_ = <IN>; print OUT $_;
    }
}
close IN;
close OUT;
