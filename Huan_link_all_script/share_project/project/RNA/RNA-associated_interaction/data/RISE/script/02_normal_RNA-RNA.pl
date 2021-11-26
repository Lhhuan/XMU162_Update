# 将../output/01_rise_human_hg19.txt 和../output/01_rise_human_position_null.txt  normal成 ../normalized/02_rise_human_hg19_normalized.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @files = ("01_rise_human_hg19.txt","01_rise_human_position_null.txt");

foreach my $file(@files){
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
        unless ($_ =~/chr1/){
            my $Chr1 =$f[0];
            my $Start1 =$f[1];
            my $End1 =$f[2];
            my $Strand1 =$f[8];
            my $Type1 = $f[14];
            my $interaction_name1 = $f[11];
            my $Gene_ID1 = $f[10];
            my $Gene_name1 = $f[11];
            my $Entrez_ID1 ="NA";
            #-----------
            my $Chr2 =$f[3];
            my $Start2 =$f[4];
            my $End2 =$f[5];
            my $Strand2 =$f[9];
            my $Type2 = $f[15];
            my $interaction_name2 = $f[13];
            my $Gene2_ID2 =$f[12];
            my $Gene_name2 =$f[13];
            my $Entrez_ID2 = "NA";
            my $Sub_interaction_type = "${Type1}_${Type2}";
            my $Interaction_type = "RNA-RNA";
            my $Genome_version = "Hg19";
            my $cell_line_or_tissue =$f[-2];
            my $pancancerNum = "NA";
            my $Source = "RISE";
            my $output = "$Chr1\t$Start1\t$End1\t$Strand1\t$Type1\t$interaction_name1\t$Gene_ID1\t$Gene_name1\t$Entrez_ID1\t$Chr2\t$Start2\t$End2\t$Strand2\t$Type2\t$interaction_name2\t$Gene2_ID2\t$Gene_name2\t$Entrez_ID2\t$Sub_interaction_type\t$Interaction_type\t$Genome_version\t$cell_line_or_tissue\t$pancancerNum\t$Source";
            print $O1 "$output\n";
        }
    }
}