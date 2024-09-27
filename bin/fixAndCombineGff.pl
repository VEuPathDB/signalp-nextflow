#!/usr/bin/perl

use strict;

use Getopt::Long;

my ($gff, $regionGff, $spVersion, $outputFile);

GetOptions("gff=s" => \$gff,
           "region_gff=s" => \$regionGff,
           "sp_version=i" => \$spVersion,
           "output_file=s" => \$outputFile
    );


open(OUT, ">$outputFile") or die "Cannot open output file $outputFile for writing: $!";

open(GFF, $gff) or die "Cannot open gff $gff file for reading: $!";
while(<GFF>) {
    chomp;
    next if /^#/;
    my @a = split(/\t/, $_);

    my ($protein) = $a[0] =~ /^(\S+)/;

    $a[0] = $protein;

    #the 9th column in gff3 is the attributes field
    $a[8] = "ID=sp${spVersion}_${protein}";

    print OUT join("\t", @a) . "\n";
}

close GFF;


if($regionGff) {

    open(REGION, $regionGff) or die "Cannot open gff $regionGff file for reading: $!";
    while(<REGION>) {
        chomp;
        next if /^#/;
        my @a = split(/\t/, $_);

        my ($protein) = $a[0] =~ /^(\S+)/;
        $a[0] = $protein;

        my $id = "sp${spVersion}_${protein}_$a[3]_$a[4]";

        #the 9th column in gff3 is the attributes field
        $a[8] = "ID=$id;Parent=sp${spVersion}_${protein}";

        print OUT join("\t", @a) . "\n";
    }
    close REGION;

}

close OUT;

1;
