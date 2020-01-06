#!/usr/bin/perl -w
# usage : JoinReferenceSegment.pl [options] ref.fa ref.mrefid.rp mrs mre

use strict;
use Getopt::Long;


##### set parameters and arguments

my $f1   = 330;     # length of reference 1
my $f2   = 330;     # length of reference 2
my $nl   = 8;       # number of Ns between paired reads for direct joining
my $join = "io";    # join method: io (inside-out), dj (direct), or amp (amplicon)

GetOptions(
    "f1=i"    => \$f1,
    "f2=i"    => \$f2,
    "nl=i"    => \$nl,
    "join=s"  => \$join,
    );

my ( $reffa, $refpos, $mrs, $mre ) = @ARGV;

my $Ns = 'N' x $nl;


### load reference positions and calculate start and end positions of amplicon

my %refas;    # amplicon start on the reference
my %refae;    # amplicon end on the reference

open IN, "<$refpos" || die "open $refpos: $!\n";

while( <IN> ) {
    my @a  = split "\t"; chomp $a[-1];
    my $id = $a[0];
    my @rp = split ',', $a[1];

    # get amplicon start position
    for( my $i=0; $i<@rp; $i++ ) {
        if( $rp[$i] > $mrs ) {
            my $d = $rp[$i] - $mrs;
            $refas{$id} = $i - $d;
            last;
        }
    }

    # get amplicon end position
    for( my $i=$#rp; $i>=0; $i-- ) {
        if( $rp[$i] < $mre ) {
            my $d = $mre - $rp[$i];
            $refae{$id} = $i + $d;
            last;
        }
    }
}
close IN;


### load reference sequence, get segment 1 and 2, and join

open IN, "<$reffa" || die "open $reffa: $!\n";
open LOG, ">$reffa.log" || die "open $reffa.log: $!\n";

while( <IN> ) {
    
    # for a new record
    if( />(\S+)/ ) {
	my $id  = $1;

	# load sequence
	my $seq;
	while( <IN> ) {
	    if( /^>/ ) {
		seek( IN, -length($_), 1 );
		last;
	    } else {
		chomp;
		$seq .= $_;
	    }
	}
	my $sl = length($seq);

	# skip if the reference does not cover either primer site, i.e., outside the amplicon
	next if !defined $refas{$id} || !defined $refae{$id};
	
	# append Ns to the left if the reference does not fully cover the forward primer
	my $si = $refas{$id};
	my $ln = "";
	if( $si < 0 ) {
	    $ln = 'N' x abs($si);
	    $si = 0;
	}

	# append Ns to the right if the reference does not fully cover the reverse primer
	my $ei = $refae{$id};
	my $rn = "";
	if( $ei > ($sl-1) ) {
	    $rn = 'N' x ($ei-$sl+1);
	    $ei = $sl - 1;
	}

	# get amplicon and length info
	my $amp = $ln . substr( $seq, $si, ($ei-$si+1) ) . $rn;
	my $al  = length($amp);
	my $lnl = length($ln);
	my $rnl = length($rn);
	my $nl  = $lnl + $rnl;

	# skip if the padded Ns constitute more than half of the amplicon
	# next if $nl > $al/2;

	# output amplicon
	if( $join eq "amp" ) {
	    print ">$id\n$amp\n";
	    print LOG "$id\t$si,$ei:$lnl,$rnl\n";
	    next;
	}

	# get reference segment
	my $l1 = $f1<=($al-$rnl) ? $f1 : $al-$rnl;
	my $l2 = $f2<=($al-$lnl) ? $f2 : $al-$lnl;
	my $seg1 = substr( $amp, 0, $l1 );
	my $seg2 = substr( $amp, -$l2 );
	my $jseg = $join eq "dj" ? $seg1.$Ns.$seg2 : $seg2.$seg1;

	# skip if either segment is full of padded Ns or the padded Ns constitute more than half of the joined segment
	# next if $lnl>=$l1 || $rnl>=$l2 || $nl > ($l1+$l2)/2;

	# output joined reference segment
	print ">$id\n$jseg\n";
	print LOG "$id\t$si,$ei:$l1,$l2:$lnl,$rnl\n";
    }
}
close IN;
