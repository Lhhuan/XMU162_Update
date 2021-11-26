# 将../output/01_rise_human_hg19.txt normal成 ../normalized/02_rise_human_hg19_normalized.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../output/$file";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../normalized/02_rise_human_hg19_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";

my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my %hash1;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless ($_ =~/start1/){
        my $Chr1 =
        my $Start1 =
        my $End1 =
        my $Strand1 =
        my $Type1 = 
        my $interaction_name1 = ;
        my $Gene_ID1 = ;
        my $Gene_name1 =;
        my $Entrez_ID1 =;
        #-----------
        my $Chr2 =;
        my $Start2 =;
        my $End2 =;
        my $Strand2 =;
        my $Type2 = ;
        my $interaction_name2 = ;
        my $Gene2_ID2 =;
        my $Gene_name2 =;
        my $Entrez_ID2 = ;
        my $Sub_interaction_type = "${Type1}_${Type2}";
        my $Interaction_type = ;
        my $Genome_version = ;
        my $cell_line_or_tissue =;
        my $pancancerNum = ;
        my $Source = ;
        my $output = "$Chr1\t$Start1\t$End1\t$Strand1\t$Type1\t$interaction_name1\t$Gene_ID1\t$Gene_name1\t$Entrez_ID1\t$Chr2\t$Start2\t$End2\t$Strand2\t$Type2\t$interaction_name2\t$Gene2_ID2\t$Gene_name2\t$Entrez_ID2\t$Sub_interaction_type\t$Interaction_type\t$Genome_version\t$cell_line_or_tissue\t$pancancerNum\t$Source";
        print $O1 "$output\n";
    }
}
