#normalize "../raw_data/ENCORI_hg19_${geneType}-RNA_RNA-RNA_all.txt", 得../normalized/08_ENCORI_hg19_RNA-RNA_normalized.txt;
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $fo1 = "../normalized/08_ENCORI_hg19_RNA-RNA_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# my $fo2 = "../output/08_ENCORI_hg19_RNA-RNA_normalized_fail.txt"; #得没有匹配的mirna文件
# open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my (%hash1,%hash2,%hash3);

my @genetypes = ("mRNA","lncRNA","pseudogene","sncRNA","miRNA");

foreach my $geneType(@genetypes){
    my $f2 = "../raw_data/ENCORI_hg19_${geneType}-RNA_RNA-RNA_all.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
    while(<$I2>) #在06中用end - start确定该文件为0-based
    {
        chomp;
        unless(/^#|geneID/){
            my @f= split/\t/;
            my $geneID1 = $f[0];
            my $geneName1 =$f[1];
            my $Type1 = $geneType;
            my $GeneID2 = $f[3];
            my $geneName2 =$f[4];
            my $Type2 = "RNA";
            my $interactionLocus =$f[10];
            my @ts = split/\|/,$interactionLocus;
            my @r1 = split/\:/,$ts[0];
            my @r2 =split/\:/,$ts[1];
            my $chr1 = $r1[0];
            my $start1 = $r1[1];
            my $end1 = $r1[2];
            my $Strand1 = $r1[3];
            my $interaction_name1 = $geneName1;
            my $chr2 = $r2[0];
            my $start2 = $r2[1];
            my $end2 = $r2[2];
            # my $ddd = $end2 -$start2 ;
            # print "$ddd\n";
            my $Strand2 = $r2[3];
            my $interaction_name2 = $geneName2;
            my $CellLine_Tissue =$f[-1];
            my $Sub_interaction_type = "${geneType}-RNA";
            my $output = "$chr1\t$start1\t$end1\t$Strand1\t$Type1\t$interaction_name1\t$geneID1\t$geneName1\tNA\t$chr2\t$start2\t$end2\t$Strand2\t$Type2\t$interaction_name2\t$GeneID2\t$geneName2\tNA\t$Sub_interaction_type\tRNA-RNA\tHg19\t$CellLine_Tissue\tNA\tstarBase";
            unless (exists $hash2{$output}){
                $hash2{$output} =1;
                print $O1 "$output\n";
            }
        }
    }
}
