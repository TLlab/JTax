#!/usr/bin/perl -w
# usage : ExtractFasta.pl fa id
# note  : This extracts the desire fasta sequence via ID.

use strict;


# extract the desired fasta sequence

open IN, "<$ARGV[0]" || die "open $ARGV[0]: $!\n";

while( <IN> ) {
    if( /^>$ARGV[1];/ ) {
	print;
	while( <IN> ) {
	    last if /^>/;
	    print;
	}
	last;
    }
}
close IN;
