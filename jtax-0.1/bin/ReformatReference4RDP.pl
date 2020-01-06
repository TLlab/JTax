#!/usr/bin/perl -w
# usage : ReformatReference4RDP.pl ref.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $out = "";

GetOptions(
    "out=s" => \$out,
    );

my $reffa = $ARGV[0];

if( !$out ) {
    my ($ref) = $reffa =~ /(.+)\.fa/;
    $out = "r$ref";
}

my @fr = qw( rootrank domain phylum class order family genus );    # fix rank


##### generate fasta and taxonomy file for retraining RDP Classifier

open OUTfa, ">$out.fa" || die "open $out.fa: $!\n";
open OUTtt, ">$out.tt" || die "open $out.tt: $!\n";

my %taxnode;

my $n = 0;
$taxnode{"Root"} = $n++;
print OUTtt "0*Root*-1*0*$fr[0]\n";

open IN, "<$reffa" || die "open $reffa: $!\n";

while( <IN> ) {
    if( /^>/ ) {
	my  ( $id, $tax ) = />(.+);tax=(.+);/;

	my @tax = split ',', $tax;
	unshift @tax, "Root";

	$tax = join( ';', @tax );
	print OUTfa ">$id\t$tax\n";

	for( my $l=1; $l<@tax; $l++) {
	    if( !$taxnode{$tax[$l]} ) {
		$taxnode{$tax[$l]} = $n++;
		print OUTtt "$taxnode{$tax[$l]}*$tax[$l]*$taxnode{$tax[$l-1]}*$l*$fr[$l]\n";
	    }
	}
    } else {
	print OUTfa $_;
    }
}
close IN;
