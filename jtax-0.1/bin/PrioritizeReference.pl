#!/usr/bin/perl -w
# usage : PrioritizeReference.pl ref.fa primer.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $maxd = 1;           # maximal number of mismatches in primer alignment
my $cpu  = 2;           # number of threads
my $maxhbn2m  = 200;    # maximal number of hanging bases to the median

GetOptions(
    "maxd=i" => \$maxd,
    "cpu=i"  => \$cpu,
    );

my ( $reffa, $primerfa ) = @ARGV;


##### prepare primer and align them to reference

my $fp = `grep -A 1 FP $primerfa | tail -1`; chomp $fp;    # forward primer
my $rp = `grep -A 1 RP $primerfa | tail -1`; chomp $rp;    # reverse primer

my $ppfa = "$primerfa.plus.fa";
open OUT, ">$ppfa" || die "open $ppfa: $!\n";
print OUT ">FP\n$fp\n";
print OUT ">RP\n" . ReverseComplement($rp) . "\n";
close OUT;

my $pn=`grep -c \\> $ppfa`; chomp $pn;

my $usearch = "usearch -threads $cpu -search_oligodb $reffa -db $ppfa -userout $reffa.$ppfa.usearch -strand plus";
$usearch .= " -userfields target+qstrand+diffs+tlor+thir+qlor+qhir+ql+qrow+trowdots+query -maxdiffs $maxd 2>/dev/null";

`$usearch`;
`rm $ppfa`;


##### load alignment result

my %reflhbn;    # number of hanging bases outside the left primer sites on each reference
my %refrhbn;    # number of hanging bases outside the right primer sites on each reference 
my %refthbn;    # total number of hanging bases outside the primer sites of each reference

open IN, "<$reffa.$ppfa.usearch" || die "open $reffa.$ppfa.usearch: $!\n";

my $lastref;
my %pbsn;              # number of binding sites of a primer
my $lhbn = 1000000;    # number of hanging bases outside the left primer
my $rhbn = 1000000;    # number of hanging bases outside the right primer

while(<IN>) {
    my @a=split("\t",$_); chomp $a[-1];

    # for a new reference
    if(!$lastref) {
	$lastref = $a[10];
        %pbsn = ();           
        $lhbn = 1000000;
        $rhbn = 1000000; 
    }

    # if the alignment is for the same reference
    if($a[10] eq $lastref) {

	# get primer binding info and count hanging bases
        $pbsn{$a[0]}++;
        my $lh = $a[5];
        my $rh = $a[7]-$a[6];
        $lhbn = $lh if $lh < $lhbn;
	$rhbn = $rh if $rh < $rhbn;
 
    # if the alignment is for a different reference
    } else {

	# select reference with a unique binding site for all primers
        my @up = grep( $pbsn{$_}==1, keys %pbsn );
        if( @up == $pn ) {
	    $reflhbn{$lastref} = $lhbn;
	    $refrhbn{$lastref} = $rhbn;
            $refthbn{$lastref} = $lhbn + $rhbn;
        }

	# reset
        seek( IN, -length($_), 1 );
        $lastref = "";
    }
}

# for the last reference
my @up = grep( $pbsn{$_}==1, keys %pbsn );
if( @up == $pn ) {
    $reflhbn{$lastref} = $lhbn;
    $refrhbn{$lastref} = $rhbn;
    $refthbn{$lastref} = $lhbn + $rhbn;
}

close IN;


# check number of hanging bases

my @hbn = sort { $b <=> $a } values %refthbn;

my $median = Median( @hbn );


# sort reference by number of hanging bases and output

my @ref = sort{ $refthbn{$b} <=> $refthbn{$a} } keys %refthbn;
my @refl = grep( $refthbn{$_} > ($median+$maxhbn2m), @ref );
my @refs = grep( $refthbn{$_} <= ($median+$maxhbn2m), @ref ); 

foreach my $r (@refs,@refl) {
    print "$r\t$refthbn{$r}\t$reflhbn{$r}\t$refrhbn{$r}\n";
}

`rm $reffa.$ppfa.usearch`;


#################### sub-routines ####################


sub ReverseComplement {
    my $s = shift;
    $s =~ tr/ACGTRYKMBVDH/TGCAYRMKVBHD/;
    $s = reverse($s);
    return $s;
}


sub Median {
    my $n = scalar( @_ );
    my $m = $n % 2 == 0 ? ($_[$n/2-1]+$_[$n/2])/2 : $_[($n-1)/2];
    return $m;
}
