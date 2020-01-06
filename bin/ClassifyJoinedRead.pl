#!/usr/bin/perl -w
# usage : ClassifyJoinedRead.pl join_ref join_read.fa/q

use strict;
use Getopt::Long;


##### set parameters and arguments

my $method = "rdp";    # taxonomy classification method
my $rdpbd  = "";       # base directory to RDP Classifier
my $ma     = 128;      # usearch taxonomy option: maximal accepts
my $mr     = 128;      # usearch taxonomy option: maximal rejects
my $idmin  = 0.84;     # usearch taxonomy option: minimal identity for usearch_global
my $pidd   = 0.5;      # usearch taxonomy option: % identity difference from maximum to accept an alignment
my $out    = "";       # output name
my $cpu    = 2;        # number of threads

GetOptions(
    "method=s"  => \$method,
    "rdpbd=s"   => \$rdpbd,
    "ma=i"      => \$ma,
    "mr=i"      => \$mr,
    "idmin=f"   => \$idmin,
    "pidd=f"    => \$pidd,
    "out=s"     => \$out,
    "cpu=i"     => \$cpu,
    );

my ( $jref, $jreadf ) = @ARGV;

my $bjref   = `basename $jref`;   chomp $bjref;
my $bjreadf = `basename $jreadf`; chomp $bjreadf;

my $f = $bjreadf =~ /a$/ ? "fa" : "fq";

$out  = "$bjref.$bjreadf.$method" if !$out;

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


##### classification in parallel

# for RDP classifier
if( $method eq "rdp" ) {

    `SplitFastx.pl $jreadf $cpu`;

    my @child=();
    
    for( my $i=1; $i<=$cpu; $i++ ) {
	
	my $pid = fork();
	if($pid) {
	    push( @child, $pid );
	    
	} elsif( $pid == 0 ) {
	    my $rdpc = "java -Xmx1g -jar $rdpbd/dist/classifier.jar classify";
	    $rdpc   .= " -t $jref.retrained/rRNAClassifier.properties -f fixrank";
	    $rdpc   .= " -o $bjref.$bjreadf.$i.rdp $jreadf.$i.$f";
	    `$rdpc`;
	    exit 0;
	    
	} else {
	    print "fork: $!\n";
	}
    }
    
    foreach (@child){
	waitpid( $_, 0 );
    }

    # combine results and remove intermediate files
    
    my $cat = "cat " . join( " ", map( "$bjref.$bjreadf.$_.$method", (1..$cpu) ) ) . " > $out";
    `$cat`;
    
    `rm $bjref.$bjreadf.*.$method`;

    my $rm = "rm " . join( " ", map( "$jreadf.$_.$f", (1..$cpu) ) );
    `$rm`;

# for other methods
} else {
    if( $method eq "sintax" ) {
	`usearch -sintax $jreadf -db $jref.udb -strand plus -threads $cpu -tabbedout $out 2> /dev/null`;

    } elsif( $method eq "ktop" ) {
	`usearch -sintax $jreadf -db $jref.udb -strand plus -threads $cpu -tabbedout $out -ktop 2> /dev/null`;
    
    } elsif( $method eq "nbc" ) {
	`usearch -nbc_tax $jreadf -db $jref.fa -strand plus -threads $cpu -tabbedout $out 2> /dev/null`;
	
    } elsif( $method eq "top" ) {
	`usearch -cons_tax $jreadf -db $jref.udb -strand plus -threads $cpu -tabbedout $out -id 0.7 -maxaccepts 3 -maxrejects 16 -top_hit_only 2> /dev/null`;

    } elsif( $method eq "btop" ) {
	if( $f eq "fa" ) {
	    `ln $jreadf $bjreadf.fa`;
	} else {
	    `Fastq2Fasta.pl $jreadf > $bjreadf.fa`;
	}
	my $blastn = "blastn -task megablast -db $jref/db -query $bjreadf.fa -num_threads $cpu ";
	$blastn .= " -max_target_seqs 1 -outfmt \"6 qseqid sseqid length pident qstart qend\" -evalue 0.01 > $bjref.$bjreadf.blast";
	`$blastn`;
	`Blast2Btop.pl $jref/db_original.fa $bjref.$bjreadf.blast > $out`;
	`rm $bjreadf.fa`;
	`rm $bjref.$bjreadf.blast`;
    }
}
