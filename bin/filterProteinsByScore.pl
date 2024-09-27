#!/usr/bin/perl

use strict;

use Getopt::Long;

use Bio::SeqIO;

my ($fasta, $predictions, $spScoreCutoff, $pctProteinsCutoff, $outputFile);

GetOptions("fasta=s" => \$fasta,
           "predictions=s" => \$predictions,
           "score_cutoff=f" => \$spScoreCutoff,
           "pct_proteins_cutoff=f" => \$pctProteinsCutoff,
           "output_file=s" => \$outputFile
    );



# NOTE:  the predictions file has headers like this which will be removed.
#        we are interested in the score and the "SP" column for filtering
# SignalP-5.0	Organism: euk	Timestamp: 20240927030521
# ID	Prediction	SP(Sec/SPI)	OTHER	CS Position

my $proteinScores = &makeProteinScoresArray($predictions);;

my $proteinCount = scalar @$proteinScores;
my $minProteins = $proteinCount * ($pctProteinsCutoff / 100);


my $keepProteins = &filterProteinsByScoreAndCount($proteinScores, $minProteins);;

&writeFilteredFasta($fasta, $keepProteins, $outputFile);


sub writeFilteredFasta {
    my ($fasta, $keepProteins) = @_;

    my $fileType = "fasta";

    my $seqio = Bio::SeqIO->new(-file => $fasta, -format => $fileType);
    my $seqout = Bio::SeqIO->new(-file => ">$outputFile", -format => $fileType);

    while (my $seq = $seqio->next_seq) {
        if($keepProteins->{$seq->id}) {
            $seqout->write_seq($seq);
        }
    }

    $seqio->close();
    $seqout->close();
}


sub filterProteinsByScoreAndCount {
    my ($proteinScores, $minProteins) = @_;

    my %keepProteins;

    # sort proteins by score desc
    my @sortedProteinScores = sort { $b->[1] <=> $a->[1] } @$proteinScores;

    my $count;
    foreach my $p (@sortedProteinScores) {
        last if($p->[1] < $spScoreCutoff && $count > $minProteins);
        $count++;
        $keepProteins{$p->[0]} = 1;
    }

    return \%keepProteins;
}


sub makeProteinScoresArray {
    my ($predictions) = @_;

    my @proteinScores;

    open(FILE, $predictions) or die "Cannot open file $predictions for reading: $!";
    while(<FILE>) {
        chomp;
        next if/^#/; #skip header lines

        my @a = split(/\t/, $_);

        push @proteinScores, [$a[0], $a[2]];
    }

    close FILE;

    return \@proteinScores;
}


1;
