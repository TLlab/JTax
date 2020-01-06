#!/usr/bin/perl -w
# usage : RearrangeReference.pl [options] ref.fa primer.fa mrefid

use strict;
use Getopt::Long;


##### set parameters and arguments

my $mrefid = "";               # main reference id
my $idmin  = 0.5;              # minimal identity of alignment between a reference and the main reference
my $l1     = 300;              # read 1 length
my $l2     = 300;              # read 2 length
my $fr     = 1;                # ratio of reference segment length to read length for inside-out rearrangement
my $trimp  = 0;                # trim primer
my $nl     = 8;                # number of Ns between paired reads for direct joining
my $join   = "io";             # joining options: io (inside-out), dj (direct joining), or amp (amplicon)
my $outfa  = "ref$join.fa";    # name of the rearranged reference
my $cpu    = 2;                # number of threads

GetOptions(
    "mrefid=s" => \$mrefid,
    "idmin=f"  => \$idmin,
    "l1=i"     => \$l1,
    "l2=i"     => \$l2,
    "fr=f"     => \$fr,
    "trimp=i"  => \$trimp,
    "nl=i"     => \$nl,
    "join=s"   => \$join,
    "outfa=s"  => \$outfa,
    "cpu=i"    => \$cpu,
    );

my ( $reffa, $primerfa ) = @ARGV;

my $breffa    = `basename $reffa`;    chomp $breffa;
my $bprimerfa = `basename $primerfa`; chomp $bprimerfa;

`ln -s $reffa` if $breffa ne $reffa;
`ln -s $primerfa` if $bprimerfa ne $primerfa;

if($join ne "io" && $join ne "dj" && $join ne "amp") {
    print STDERR "Error: only io (inside-out), dj (direct joining), or amp (amplicon) is accpeted as a joining option!\n";
    exit 0;
}

my $f1 = $join eq "d" ? $l1 : int( $l1 * $fr );
my $f2 = $join eq "d" ? $l2 : int( $l2 * $fr );


##### get main reference and align primer to the main reference

if( !$mrefid ) {
    `PrioritizeReference.pl $breffa $bprimerfa > $breffa.$bprimerfa.hbn`;
    my $maxhbn = `head -1 $breffa.$bprimerfa.hbn`;
    ($mrefid) = $maxhbn =~ /(.+);tax/;
    #($mrefid) = $maxhbn =~ /^(\S+)/;
    `rm $breffa.$bprimerfa.hbn`;
}

`ExtractFasta.pl $breffa $mrefid > $mrefid.fa`;    # extract main reference sequence via ID   

# rearrange primer to be on the plus strand of reference

my $fp = `grep -A 1 FP $bprimerfa | tail -1`; chomp $fp;
my $rp = `grep -A 1 RP $bprimerfa | tail -1`; chomp $rp;

open OUT, ">$bprimerfa.plus.fa" || die "open $bprimerfa.plus.fa: $!\n";
print OUT ">FP\n" . "$fp\n";
print OUT ">RP\n" . ReverseComplement($rp) . "\n";
close OUT;

# do alignment (possibility of more than one primer site?)

my $usearch = "usearch -search_oligodb $mrefid.fa -db $bprimerfa.plus.fa -threads $cpu -userout $mrefid.fa.$bprimerfa.usearch -strand plus";
$usearch .= " -userfields target+qstrand+diffs+tlor+thir+tl+qlor+qhir+ql+qrow+trowdots+query -maxdiffs 1 2> /dev/null";
`$usearch`;

`rm $bprimerfa.plus.fa`;


##### get positions of the forward and reverse primers on the main reference

my $mrs;    # main reference start at the forward primer
my $mre;    # main reference end at the reverse primer

open IN, "<$mrefid.fa.$bprimerfa.usearch" || die "open $mrefid.fa.$bprimerfa.usearch: $!\n";

while( <IN> ) {
    my @a = split "\t";

    if($a[0] eq "FP") {
	$mrs = $trimp==0 ? $a[6]+1 : $a[7]+2;

    } elsif($a[0] eq "RP") {
	$mre = $trimp==0 ? $a[7]+1 : $a[6];
    }
}
close IN;

`rm $mrefid.fa.$bprimerfa.usearch`;

if( !$mrs || !$mre ) {
    print STDERR "Error: forward or reverse primer not found in the main reference!\n";
    exit 0;
}


##### split reference sequences for parallel processing

`SplitFastx.pl $breffa $cpu`;

my @child = ();

for(my $i=1; $i<=$cpu; $i++) {

    my $pid = fork();
    if($pid) {
        push( @child, $pid );
    } elsif( $pid == 0 ) {

        # map base position to the main reference
	`MapBasePosition2MainReference.pl -idmin $idmin $breffa.$i.fa $mrefid.fa > $breffa.$mrefid.$i.rp`;

        # process sequence into inside-out reference
	`JoinReferenceSegment.pl -f1 $f1 -f2 $f2 -nl $nl -j $join $breffa.$i.fa $breffa.$mrefid.$i.rp $mrs $mre > $breffa.$mrefid.$i.$join.fa`;

        exit 0;

    } else {
        print "fork: $!\n";
    }
}

foreach (@child) {
    waitpid $_, 0;
}


##### combine joined references and remove intermediate files

my $cat = "cat " . join( " ", map( "$breffa.$mrefid.$_.$join.fa", (1..$cpu) ) ) . " > $outfa";
`$cat`;

my $rm = "rm " . join( " ", map( "$breffa.$mrefid.$_.$join.fa", (1..$cpu) ) );
`$rm`;

$cat = "cat " . join( " ", map( "$breffa.$_.fa.log", (1..$cpu) ) ) . " > $outfa.log";
`$cat`;

$rm = "rm " . join( " ", map( "$breffa.$_.fa.log", (1..$cpu) ) );
`$rm`;

$rm = "rm " . join( " ", map( "$breffa.$_.fa", (1..$cpu) ) );
`$rm`;

$rm = "rm " . join( " ", map( "$breffa.$mrefid.$_.rp", (1..$cpu) ) );
`$rm`;


`rm $mrefid.fa`;

`rm $breffa` if $breffa ne $reffa;
`rm $bprimerfa` if $bprimerfa ne $primerfa;



######################################################################


sub ReverseComplement {
    my $s = shift;
    $s =~ tr/acgtrykmbvdhACGTRYKMBVDH/tgcayrmkvbhdTGCAYRMKVBHD/;
    $s = reverse($s);
    return $s;
}
