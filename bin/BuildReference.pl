#!/usr/bin/perl -w
# usage : BuildReference.pl [option] ref.fa primer.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $mrefid = "";            # main reference id
my $idmin  = 0.5;           # minimal identity of alignment between a reference and the main reference
my $l1     = 300;           # read 1 length
my $l2     = 300;           # read 2 length
my $fr     = 1;             # ratio of reference segment length to read length for inside-out rearrangement
my $trimp  = 0;             # trim primer
my $nl     = 8;             # number of Ns between paired reads for direct joining
my $join   = "io";          # joining options: io (inside-out), dj (direct joining), or amp (amplicon)
my $out    = "ref$join";    # base name of the rearranged output
my $method = "rdp";         # supported methods: rdp, nbc, sintax, ktop, top, btop
my $rdpbd  = "";            # base directory to RDP Classifier
my $cpu    = 2;             # number of threads

GetOptions(
    "mrefid=s"  => \$mrefid,
    "idmin=f"   => \$idmin,
    "l1=i"      => \$l1,
    "l2=i"      => \$l2,
    "fr=f"      => \$fr,
    "trimp=i"   => \$trimp,
    "nl=i"      => \$nl,
    "join=s"    => \$join,
    "out=s"     => \$out,
    "method=s"  => \$method,
    "rdpbd=s"   => \$rdpbd,
    "cpu=i"     => \$cpu,
    );

my ( $reffa, $primerfa ) = @ARGV;

# uncomment the following line and set the absolute path to avoid the long input
$rdpbd = "/home/tliu/tool/rdpclassifier/rdp_classifier_2.12";    # base directory to RDP classifier


##### rearrange reference sequences for the primer pair

my $raref = "RearrangeReference.pl -i $idmin -l1 $l1 -l2 $l2 -fr $fr -t $trimp -nl $nl -j $join -o $out.fa";
$raref .= " -mrefid $mrefid" if $mrefid;
$raref .= " $reffa $primerfa";
`$raref`;


##### create database using the rearranged reference

my $cd = "CreateDatabase.pl -m $method -rdpbd $rdpbd -o $out $out.fa";
`$cd`;

print "$cd\n";
