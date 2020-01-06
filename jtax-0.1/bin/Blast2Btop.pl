#!/usr/bin/perl -w
# usage : BtopRaw2Pred.pl train.fa btop.raw

use strict;


# load taxonomy

my %acctax;

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while( <IN> ) {
    if( /^>(.+);tax=(.+);/ ) {
	$acctax{$1}=$2;
    }
}
close IN;


# load btop results

my %ridq;

open IN, "<$ARGV[1]" || die "open $ARGV[1]: $!\n";

while( <IN> ) {
    my @a = split "\t"; chomp $a[-1];
    if(!$ridq{$a[0]}) {
	my $acc = lc( $a[1] );
	print "$a[0]\t$acctax{$acc}\t".join("\t",@a[2..$#a])."\n";
	$ridq{$a[0]}=1;
    }
}
close IN;
