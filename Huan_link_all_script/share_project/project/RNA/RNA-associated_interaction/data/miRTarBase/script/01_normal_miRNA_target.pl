# 将../raw_data/miRTarBase_MTI.txt normal成 ../normalized/01_miRTarBase_normalized.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../raw_data/miRTarBase_MTI.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../normalized/01_miRTarBase_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";

my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my %hash1;
while(<$I1>)
{
    chomp;
    unless ($_ =~/miRTarBase/){
        my @f= split/\t/;
        my $miRNA =$f[1];
        my $Species_miRNA = $f[2];
        my $Target_Gene =$f[3];
        my $Target_Gene_Entrez_ID = $f[4];
        my $Species_Target_Gene = $f[5];
        if($Species_miRNA =~/Homo sapiens/ && $Species_Target_Gene=~/Homo sapiens/){
            my $Chr1 = "NA";
            my $Start1 ="NA";
            my $End1 ="NA";
            my $Strand1 ="NA";
            my $Type1 = "miRNA";
            my $interaction_name1 = $miRNA;
            my $Gene_ID1 = "NA";
            my $Gene_name1 ="NA";
            my $Entrez_ID1 ="NA";
            #-----------
            my $Chr2 ="NA";
            my $Start2 ="NA";
            my $End2 ="NA";
            my $Strand2 ="NA";
            my $Type2 ="RNA" ;
            my $interaction_name2 =$Target_Gene ;
            my $Gene2_ID2 ="NA";
            my $Gene_name2 =$Target_Gene;
            my $Entrez_ID2 =$Target_Gene_Entrez_ID ;
            my $Sub_interaction_type = "${Type1}_${Type2}";
            my $Interaction_type = "RNA-RNA";
            my $Genome_version = "NA";
            my $cell_line_or_tissue ="NA";
            my $pancancerNum = "NA";
            my $Source ="miRTarBase" ;
            my $output = "$Chr1\t$Start1\t$End1\t$Strand1\t$Type1\t$interaction_name1\t$Gene_ID1\t$Gene_name1\t$Entrez_ID1\t$Chr2\t$Start2\t$End2\t$Strand2\t$Type2\t$interaction_name2\t$Gene2_ID2\t$Gene_name2\t$Entrez_ID2\t$Sub_interaction_type\t$Interaction_type\t$Genome_version\t$cell_line_or_tissue\t$pancancerNum\t$Source";
            print $O1 "$output\n";
        }
    }
}
