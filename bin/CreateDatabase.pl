#!/usr/bin/perl -w
# usage : CreateDatabase.pl [options] ref.fa

use strict;
use Getopt::Long;


##### set parameters and arguments

my $method = "rdp";      # taxonomy classification method
my $rdpbd  = "";         # base directory to RDP Classifier
my $out    = "refio";    # base name of the output database

GetOptions(
    "method=s" => \$method,
    "rdpbd=s"  => \$rdpbd,
    "out=s"    => \$out,
    );

my $reffa = $ARGV[0];

if( $method ne "rdp" && $method ne "nbc" && $method ne "sintax" &&
    $method ne "ktop" && $method ne "top" && $method ne "btop" ) {
    print STDERR "Error: $method not supported!\n";
    exit 0;
}

# uncomment the following line and set the absolute path to avoid the long input
$rdpbd = "/home/tliu/tool/rdpclassifier/rdp_classifier_2.12";    # base directory to RDP classifier

if( $method eq "rdp" ) {
    if( !$rdpbd ) {
	print STDERR "Error: base directory to RDP Classifier not set\n";
	exit 0;

    } elsif( !-e "$rdpbd/dist/classifier.jar" ) {
        print STDERR "Error: RDP Classifier not found\n";
        exit 0;
    }
}


##### create reference database

# for RDP Classifier
if($method eq "rdp") {

    # reformat reference data for retraining RDP Classifier
    my $breffa = `basename $reffa`;
    my ($ref) = $breffa =~ /(.+)\.fa/;
    my $rref = "r$ref";
    `ReformatReference4RDP.pl -o $rref $reffa`;
    
    # retrain RDP Classifier
    my $retrain = "java -Xmx2g -jar $rdpbd/dist/classifier.jar train";
    $retrain .= " -o $out.retrained -s $rref.fa -t $rref.tt";
    `$retrain`;
    `cp $rdpbd/samplefiles/rRNAClassifier.properties $out.retrained/`;

    # remove intermediate files, i.e., reformatted reference
    `rm $rref.fa`;
    `rm $rref.tt`;
    
# for BLAST
} elsif( $method eq "btop" ) {

    `mkdir $out`;
    `cp $reffa $out/db_original.fa`;
    
    my $sed = "sed \"-es/;.*//\" < $reffa > $out/db.fa";
    `$sed`;
    
    `makeblastdb -in $out/db.fa -dbtype nucl -parse_seqids -out $out/db`;
    
# for methods involving USEARCH except for nbc
} else {
    if($method ne "nbc") {
	`usearch -makeudb_usearch $reffa -output $out.udb 2> /dev/null`;
    }
}
