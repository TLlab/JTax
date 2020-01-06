#!/usr/bin/perl -w
# usage : jtax.pl [options] ref.fa primer.fa r1.fq r2.fq
# or    : jtax.pl [options] ref.fa primer.fa read.fa/q

use strict;
use Getopt::Long;


#################### set parameters and arguments ####################

my $refjoin  = "";       # rearranged reference database
my $method   = "rdp";    # supported methods: rdp, nbc, sintax, ktop, top, btop
my $mrefid   = "";       # main reference id
my $idmin    = 0.5;      # minimal identity of alignment between a reference and the main reference
my $l1       = 300;      # desired read 1 length
my $l2       = 300;      # desired read 2 length
my $fr       = 1;        # ratio of reference segment length to read length for inside-out rearrangement
my $trimp    = 0;        # trim primer
my $join     = "io";     # joining options: io (inside-out) or dj (direct joining)
my $minovlen = 16;       # minimal overlap length
my $pctid    = 75;       # minimal identity of overlap alignment
my $nl       = 8;        # number of Ns between paired reads for direct joining
my $out      = "";       # name of the output folder
my $cpu      = 2;        # number of threads

GetOptions(
    "refjoin=s"  => \$refjoin,
    "method=s"   => \$method,
    "mrefid=s"   => \$mrefid,
    "idmin=f"    => \$idmin,
    "l1=i"       => \$l1,
    "l2=i"       => \$l2,
    "fr=f"       => \$fr,
    "trimp=i"    => \$trimp,
    "nl=i"       => \$nl,
    "join=s"     => \$join,
    "minovlen=i" => \$minovlen,
    "pctid=f"    => \$pctid,
    "out=s"      => \$out,
    "cpu=i"      => \$cpu,
    );

if( @ARGV<3 || @ARGV>4 ) {
    Usage();
    exit 0;
}

my ( $reffa, $primerfa, $r1f ) = @ARGV;
my $r2fq = $ARGV[3] ? $ARGV[3] : "";

# basenames of files
my $breffa    = `basename $reffa`;         chomp $breffa;
my $bprimerfa = `basename $primerfa`;      chomp $bprimerfa;
my $br1f      = `basename $r1f`;           chomp $br1f;
my $br2fq     = "";
if( $r2fq ) {
    $br2fq = `basename $r2fq`; chomp $br2fq;
}

my ($ref) = $breffa =~ /(.+)\.fa/;
my ($read) = $br1f =~ /(.+)\.f/;
while( $read =~ /[_1]$/ ) { chop $read; }

# uncomment the following line and set the absolute path to avoid the long input if using RDP Classifier
my $rdpbd = "/home/tliu/tool/rdpclassifier/rdp_classifier_2.12";    # base directory to RDP Classifier


#################### verify parameters and arguments ####################

# check method

if( $method ne "rdp" && $method ne "nbc" && $method ne "sintax" &&
    $method ne "ktop" && $method ne "top" && $method ne "btop" ) {
    print STDERR "Error: $method not supported!\n";
    exit 0;

} elsif( $method eq "rdp" ) {
    if( !$rdpbd ) {
	print STDERR "Error: base directory to RDP Classifier not set in jtax.pl\n";
	exit 0;

    } elsif( !-e "$rdpbd/dist/classifier.jar" ) {
	print STDERR "Error: RDP Classifier not found\n";
	exit 0;
    }
    
} elsif( $method eq "btop" ) {
    my $blastnh = `blastn -h`;
    if( $blastnh !~ /USAGE/ ) {
	print STDERR "Error: BLAST not available\n";
	exit 0;
    }

} else {
    my $usearch = `usearch`;
    if( $usearch !~ /Edgar/ ) {
	print STDERR "Error: USEARCH not available\n";
	exit 0;
    }
}

# check reference

if( $refjoin ) {
    my $refq = 1;

    if( $method eq "rdp" ) {
	$refq = 0 if !-d "$refjoin.retrained";

    } elsif( $method eq "btop" ) {
    	$refq = 0 if !-e "$refjoin.nhr";

    } elsif( $method eq "nbc" ) {
	$refq = 0 if !-e "$refjoin.fa";

    } else {
	$refq = 0 if !-e "$refjoin.udb";
    }

    if( $refq == 0 ) {
	print STDERR "Error: reference database $refjoin not available\n";
	exit 0;
    }

} else {
    # check ID line format in the first line of the reference file
    if( !-e $reffa ) {
	print STDERR "Error: $reffa not available\n";
	exit 0;
	
    } else {
	my $fl = `head -1 $reffa`;
	if( $fl !~ /^>(.+);tax=.+;/ ) {
	    print STDERR "Error: ID line of $reffa should be >accession;tax=ranked_taxonomy;\n";
	    exit 0;
	    
	} elsif( $mrefid ) {
	    my $g = `grep $mrefid $reffa`;
	    if( !$g ) {
		print STDERR "Error: $mrefid not found in $breffa\n";
		exit 0;
	    }
	}
    }
}

# check primer file

if( !-e $primerfa ) {
    print STDERR "Error: $primerfa not available!\n";
    exit 0;

} else {
    my $gfp = `grep FP $primerfa`; chomp $gfp;
    my $grp = `grep RP $primerfa`; chomp $grp;
    if( !$gfp && !$grp ) {
	print STDERR "Error: forward and reverse primer IDs should contain FP and RP respectively\n";
	exit 0;
    }
}

# check read file

if( $r2fq ) {
    if( !-e $r1f || !-e $r2fq ) {
	print STDERR "Error: $r1f or $r2fq not available!\n";
	exit 0;
	
    } else {
	if( $r1f !~ /q$/ || $r2fq !~ /q$/ ) {
	    print STDERR "Error: $r1f and $r2fq should be in FASTQ format\n";
	    exit 0;
	}
    }
} else {
    if( $r1f !~ /a$/ && $r1f !~ /q$/ ) {
	print STDERR "Error: $r1f should be in FASTA/Q format\n";
	exit 0;
    }
}

# check other parameters

my $pq = 1;

$pq = 0 if $idmin > 1;
$pq = 0 if $l1 < 100 || $l2 < 100;
$pq = 0 if $trimp != 0 && $trimp != 1;
$pq = 0 if $join ne "io" || $join ne "dj";
$pq = 0 if $minovlen < 10;
$pq = 0 if $pctid < 50 || $pctid > 100;

Usage() if $pq;

# check output folder

if( !$out ) {
    $out = $refjoin ? $refjoin : "$ref$join";
    $out .= "\_$read\_$method";
}

if( -d $out ) {
    print "STDERR: output folder $out exists\n";
    exit 0;

} else {
    `mkdir $out`;
    chdir($out);

    `ln -s ../$reffa` if !$refjoin;
    `ln -s ../$primerfa` if !$refjoin || $r2fq;
    `ln -s ../$r1f`;
    `ln -s ../$r2fq` if $r2fq;
}


#################### run JTax pipeline ####################

# build reference

if( !$refjoin ) {
    my $br = "BuildReference.pl";
    $br .= " -mrefid $mrefid" if $mrefid;
    $br .= " -idmin $idmin -l1 $l1 -l2 $l2 -trimp $trimp -nl $nl -join $join -out $ref$join";
    $br .= " -method $method -rdpbd $rdpbd -cpu $cpu $breffa $bprimerfa";
    
    `$br`;
}

# join paired reads

if( $r2fq ) {
    my $jpr = "JoinPairedRead.pl -l1 $l1 -l2 $l2 -minovlen $minovlen -pctid $pctid -trimp $trimp";
    $jpr .= " -join $join -nl $nl -out $read.$join -cpu $cpu $br1f $br2fq $bprimerfa";
    
    `$jpr`;
}

# classify joined reads

$out =~ s/\_$method/\.$method/;

my $cjr = "ClassifyJoinedRead.pl -method $method -rdpbd $rdpbd -cpu $cpu";
$cjr .= $refjoin ? " -out $out ../$refjoin" : " -out $out $ref$join";
$cjr .= $r2fq ? " $read.$join.fq" : " $br1f";

`$cjr`;

chdir("..");


#################### sub-routines ####################


sub Usage {
    print "usage  : jax.pl [options] ref.fa primer.fa read1.fq read2.fq\n";
    print "or     : jax.pl [options] ref.fa primer.fa read.fa/q\n\n";
    print "option : -method      <str> taxonomy classification method                                              [rdp*, nbc, sintax, ktop, top, btop]\n";
    print "       : -mrefid      <str> main reference ID                                                           []\n";
    print "       : -idmin       <flt> minimal identity of alignment between a reference and the main reference    [0.5; <=1]\n";
    print "       : -l1          <int> desired read 1 length                                                       [300; >=100]\n";
    print "       : -l2          <int> desired read 2 length                                                       [300; >=100]\n";
    print "       : -trimp       <int> trim primer or not                                                          [1*,0]\n";
    print "       : -join        <str> joining method                                                              [io*,dj]\n";    
    print "       : -minovlen    <int> minimal overlap length                                                      [16; >=10]\n";
    print "       : -pctid       <flt> minimal identity of overlap alignment                                       [0.75; <=1]\n";
    print "       : -nl          <int> length of padded Ns for direct joining                                      [8]\n";
    print "       : -out         <str> name of output folder                                                       []\n";
    print "       : -cpu         <int> number of threads                                                           [2]\n";
    print "       : -help              help\n\n";
}
