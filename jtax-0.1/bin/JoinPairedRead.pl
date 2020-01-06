#!/usr/bin/perl -w
# usage : JoinPairedRead.pl [options] r1.fq r2.fq primer.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $l1       = 300;     # trim read 1 to length if longer
my $l2       = 300;     # trim read 2 to length if longer
my $minovlen = 16;      # minimal overlap length
my $pctid    = 75;      # minimal identity of overlap alignment
my $trimp    = 0;       # trim primer
my $ecq      = 0;       # error correction via overlap
my $nl       = 8;       # number of Ns between paired reads for direct joining
my $join     = "io";    # joining options: io (inside-out) or d (direct joining)
my $out      = "out";   # output file prefix
my $cpu      = 2;       # threads

GetOptions(
    "l1=i"       => \$l1,
    "l2=i"       => \$l2,
    "minovlen=i" => \$minovlen,
    "pctid=f"    => \$pctid,
    "trimp=i"    => \$trimp,
    "ecq=i"      => \$ecq,
    "nl=i"       => \$nl,
    "join=s"     => \$join,
    "out=s"      => \$out,
    "cpu=i"      => \$cpu,
    );

my ( $r1fq, $r2fq, $primerfa ) = @ARGV;

my $br1fq = `basename $r1fq`; chomp $br1fq;
my $br2fq = `basename $r2fq`; chomp $br2fq;
`ln -s $r1fq $br1fq` if $r1fq ne $br1fq;
`ln -s $r2fq $br2fq` if $r2fq ne $br2fq;

my $minmergelen = $l1;                              # minimal length for merged read
my $maxdiffs    = $minmergelen * ( 1 - $pctid );    # maximal no. of mismatches in overlap

my $Ns = 'N' x $nl;
my $Is = 'I' x $nl;


##### trim primer from paired-end reads

if( $trimp != 0 ) {

    # prepare primer file
    `grep -A 1 FP $primerfa > fp.fa`;
    `grep -A 1 RP $primerfa > rp.fa`;
    
    # align primer to read 1
    my $usearch = "usearch -search_oligodb $br1fq -db fp.fa -threads $cpu -userout $br1fq.fp.usearch -strand plus";
    $usearch .= " -userfields target+query+qstrand+diffs+tlor+thir+tl+qlor+qhir+ql -maxdiffs 3 2> /dev/null";
    `$usearch`;
    
    # get forward primer position
    my %qfpe;
    open IN, "<$br1fq.fp.usearch" || die "open $br1fq.fp.usearch: $!\n";
    while( <IN> ) {
	my @a = split "\t";
	next if $a[0]!~/FP/; 
	$qfpe{$a[1]} = $a[8]+1;
    }
    close IN;
    
    # repeat for read 2
    $usearch = "usearch -search_oligodb $br2fq -db rp.fa -threads $cpu -userout $br2fq.rp.usearch -strand plus";
    $usearch .= " -userfields target+query+qstrand+diffs+tlor+thir+tl+qlor+qhir+ql -maxdiffs 3 2> /dev/null";
    `$usearch`;
    
    my %qrpe;
    open IN, "<$br2fq.rp.usearch" || die "open $br2fq.rp.usearch: $!\n";
    while( <IN> ) {
	my @a = split "\t";
	next if $a[0]!~/RP/;
	$qrpe{$a[1]} = $a[8]+1;
    }
    close IN;

    # remove intermediate files
    `rm fp.fa`;
    `rm rp.fa`;
    `rm $br1fq.fp.usearch`;
    `rm $br2fq.rp.usearch`;
    
    # trim primers
    open IN1, "<$br1fq" || die "open $br1fq: $!\n";
    open IN2, "<$br2fq" || die "open $br2fq: $!\n";
    open OUT1, ">$br1fq.pt.fq" || die "open $br1fq.pt.fq: $!\n";
    open OUT2, ">$br2fq.pt.fq" || die "open $br2fq.pt.fq: $!\n";
    
    while( <IN1> ) {
	my $id1 = $_;
	my $s1  = <IN1>; <IN1>;
	my $q1  = <IN1>;
	
	my $id2 = <IN2>;
	my $s2  = <IN2>; <IN2>;
	my $q2  = <IN2>;
	
	# if forward and reverse primers are found in read 1 and 2 respectively
	my ($q) = $id1 =~ /^\@(\S+)/;
	
	if( $qfpe{$q} && $qrpe{$q} ) {
	    $s1 = substr( $s1, $qfpe{$q} );
	    $q1 = substr( $q1, $qfpe{$q} );
	    $s2 = substr( $s2, $qrpe{$q} );
	    $q2 = substr( $q2, $qrpe{$q} );
	    
	    print OUT1 "$id1$s1\+\n$q1";
	    print OUT2 "$id2$s2\+\n$q2";
	}
    }
    close IN1;
    close IN2;
    close OUT1;
    close OUT2;
    
} else {
    `ln -s $br1fq $br1fq.pt.fq`;
    `ln -s $br2fq $br2fq.pt.fq`;
}


##### correct errors via overlapping PE reads

if( $ecq ) {
    
    # first merge PE reads using usearch
    
    `usearch -fastq_mergepairs $br1fq.pt.fq -reverse $br2fq.pt.fq -fastqout $br1fq.$br2fq.ptmerged.fq -fastq_minovlen $minovlen -fastq_pctid $pctid -fastq_maxdiffs $maxdiffs 2> /dev/null`;
    
    # load merged results
    
    my %idms;
    my %idmq;
    
    open IN, "<$br1fq.$br2fq.ptmerged.fq" || die "open $br1fq.$br2fq.ptmerged.fq: $!\n";
    while( <IN> ) {
	my ($id) = /^\@(\S+)/;
	$idms{$id} = <IN>; chomp $idms{$id}; <IN>;
	$idmq{$id} = <IN>; chomp $idmq{$id};
    }
    close IN;
    
    # output corrected data if available
    
    open IN1, "<$br1fq.pt.fq" || die "open $br1fq.pt.fq: $!\n";
    open IN2, "<$br2fq.pt.fq" || die "open $br2fq.pt.fq: $!\n";
    open OUT1, ">$br1fq.ptec.fq" || die "open $br1fq.ptec.fq: $!\n";
    open OUT2, ">$br2fq.ptec.fq" || die "open $br2fq.ptec.fq: $!\n";
    
    while( <IN1> ) {
	my ($id) = /^\@(\S+)/;
	
	# if the pair can be merged, extract paired reads from the merged read
	if( $idms{$id} ) {
	    my $s1 = substr( $idms{$id}, 0, $l1 );
	    my $q1 = substr( $idmq{$id}, 0, $l1 );
	    my $s2 = ReverseComplement( substr( $idms{$id}, -$l2 ) );
	    my $q2 = reverse( substr( $idmq{$id}, -$l2 ) );
	    
	    print OUT1 "\@$id\n$s1\n+\n$q1\n";
	    print OUT2 "\@$id\n$s2\n+\n$q2\n";
	    
	    <IN1>; <IN1>; <IN1>;
	    <IN2>; <IN2>; <IN2>; <IN2>;
	    
	    # if the pair cannot be merged, output the original paired reads
	} else {
	    print OUT1 $_;
	    $_ = <IN1>; print OUT1 $_;
	    $_ = <IN1>; print OUT1 $_;
	    $_ = <IN1>; print OUT1 $_;
	    
	    $_ = <IN2>; print OUT2 $_;
	    $_ = <IN2>; print OUT2 $_;
	    $_ = <IN2>; print OUT2 $_;
	    $_ = <IN2>; print OUT2 $_;
	}
    }
    close IN1;
    close IN2;
    close OUT1;
    close OUT2;
    
    `rm $br1fq.$br2fq.ptmerged.fq`;

} else {
    `ln -s $br1fq.pt.fq $br1fq.ptec.fq`;
    `ln -s $br2fq.pt.fq $br2fq.ptec.fq`;
}


##### load data and join paired reads

# for inside-out joining
if( $join eq "io" ) {
    open IN1, "<$br1fq.ptec.fq" || die "open $br1fq.ptec.fq: $!\n";
    open IN2, "<$br2fq.ptec.fq" || die "open $br2fq.ptec.fq: $!\n";
    open OUT, ">$out.fq" || die "open $out.fq: $!\n";
	
    while( <IN1> ) {
	
	# load paired reads 
	my ($sid) = /^@(\S+)/;
	$sid =~ s/\/1$//;
	
	my $s1=<IN1>; chomp $s1; <IN1>;
	my $q1=<IN1>; chomp $q1;
	
	<IN2>;
	my $s2=<IN2>; chomp $s2; <IN2>;
	my $q2=<IN2>; chomp $q2;
	
	# trim to the desired length
	my $r1l = $l1<length($s1) ? $l1 : length($s1); 
	my $r2l = $l2<length($s2) ? $l2 : length($s2);
	$s1 = substr( $s1, 0, $r1l );
	$s2 = substr( $s2, 0, $r2l );
	$q1 = substr( $q1, 0, $r1l );
	$q2 = substr( $q2, 0, $r2l );
	
	# generate joined fastq
	my $seq  = ReverseComplement($s2) . $s1;
	my $qual = reverse($q2) . $q1;
	print OUT "\@$sid\n$seq\n+\n$qual\n";
    }
    close IN1;
    close IN2;

# for direct joining
} else {
    `usearch -fastq_join $br1fq.ptec.fq -reverse $br2fq.ptec.fq -fastqout $out.fq 2> /dev/null`;
}


##### remove intermediate files

`rm $br1fq.pt.fq`;
`rm $br2fq.pt.fq`;
`rm $br1fq.ptec.fq`;
`rm $br2fq.ptec.fq`;

`rm $br1fq` if $r1fq ne $br1fq;
`rm $br2fq` if $r2fq ne $br2fq;



#################### sub-routine ####################


sub ReverseComplement {
    my $s = shift;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    $s = reverse($s);
    return $s;
}
