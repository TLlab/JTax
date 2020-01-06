#!/usr/bin/perl -w 
# usage : MapBasePosition2MainReference.pl ref.fa mrefid.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $idmin = 0.5;    # minimal identity of alignment between a reference and the main reference

GetOptions(
    "idmin=f" => \$idmin,
    );

my ( $reffa, $mrefidfa ) = @ARGV;


##### align sequences to the main reference using USEARCH

my $usearch = "usearch -threads 1 -usearch_global $reffa -db $mrefidfa -strand plus -id $idmin -maxaccepts 0 -maxrejects 0 -alnout $reffa.$mrefidfa.usearch 2> /dev/null";
`$usearch`;


##### map base position to the main reference

open IN, "<$reffa.$mrefidfa.usearch" || die "open $reffa.$mrefidfa.usearch: $!\n";

while( <IN> ) {
    if( /^Query >(\S+)/ ) {
        my $q     = $1;
        my $qaseq = "";    # aligned query sequence
        my $taseq = "";    # aligned target sequence (i.e., the main reference)

	# get query length
	<IN>; <IN>; <IN>;
	$_ = <IN>;
	my ($ql) = /Query\s+(\d+)nt/;

        # get first line of alignment
	# qs/qe : query start/end
	# ts/te : target start/end
        <IN>; <IN>;
        $_ = <IN>;
        my ( $qs, $qaln, $qe ) = /^Qry\s+(\d+) \+ (\S+) (\d+)$/;
        $qaseq .= $qaln;
        <IN>;
        $_ = <IN>;
        my ( $ts, $taln, $te ) = /^Tgt\s+(\d+) \+ (\S+) (\d+)$/;
        $taseq .= $taln;
        <IN>;

        # get rest lines of alignment
        while( <IN> ) {
            last if /cols/;

            ( $qaln, $qe ) = /^Qry\s+\d+ \+ (\S+) (\d+)$/;
            $qaseq .= $qaln;
            <IN>;
            $_ = <IN>;
            ( $taln, $te ) = /^Tgt\s+\d+ \+ (\S+) (\d+)$/;
            $taseq .= $taln;
            <IN>;
        }

        # calculate target position for query
        my @qtp;

	# add position before alignment start
	for( my $i=1; $i<$qs; $i++ ) {
	    push( @qtp, $ts+($i-$qs) );
	}

	# add aligned position
        my $qp = $qs - 1;
        my $tp = $ts - 1;
        my @qab = split "", $qaseq;
        my @tab = split "", $taseq;
	
        for( my $i=0; $i<@qab; $i++ ) {
            if( $tab[$i] ne '-' ) { $tp++; }
            if( $qab[$i] ne '-' ) {
                $qp++;
                push( @qtp, $tp );
            }
        }

	# add position after alignment end
	for( my $i=$qe+1; $i<=$ql; $i++ ) {
	    push( @qtp, $te+($i-$qe) );
	}
		
        # output base position on the main reference for each query
        print "$q\t" . join(",",@qtp) . "\n";
    }
}
close IN;

# remove usearch result

`rm $reffa.$mrefidfa.usearch`;
