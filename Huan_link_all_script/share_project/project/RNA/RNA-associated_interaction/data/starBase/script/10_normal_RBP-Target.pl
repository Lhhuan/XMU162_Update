#normalize "../raw_data/ENCORI_hg19_RBP-${geneType}_RBP-Target.txt", 得../normalized/10_ENCORI_hg19_RBP-Target_normalized.txt;
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $fo1 = "../normalized/10_ENCORI_hg19_RBP-Target_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# my $fo2 = "../output/08_ENCORI_hg19_RNA-RNA_normalized_fail.txt"; #得没有匹配的mirna文件
# open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my (%hash1,%hash2,%hash3);

my @genetypes = ("mRNA","lncRNA","pseudogene","circRNA","sncRNA");

foreach my $geneType(@genetypes){
    my $f2 = "../raw_data/ENCORI_hg19_RBP-${geneType}_RBP-Target.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
    while(<$I2>) #06中用end - start确定该文件为0-based
    {
        chomp;
        unless(/^#|RBP/){
            my @f= split/\t/;
            my $protein = $f[0];
            my $interaction_name2 = $protein;

            my $Type1 = $geneType;
            #--------------
            my $GeneID2 = "NA";
            my $geneName2 ="NA";
            my $Type2 = "RBP";
            #--------------
            my $chr2 = "NA";
            my $start2 = "NA";
            my $end2 = "NA";
            my $Strand2 ="NA";
            #-------------------
            my $geneID1 = $f[1];
            my $geneName1=$f[2];
            my $interaction_name1 = $geneName1;
            my $chr1 = $f[8];
            my $start1 = $f[11];
            my $end1 = $f[12];
            my $Strand1 = $f[13];
            my $CellLine_Tissue ="NA";
            my $Sub_interaction_type = "RBP-${geneType}";
            if ($geneType =~/circRNA/){
                my $pancancerNum = "NA";
                my $output = "$chr1\t$start1\t$end1\t$Strand1\t$Type1\t$interaction_name1\t$geneID1\t$geneName1\tNA\t$chr2\t$start2\t$end2\t$Strand2\t$Type2\t$interaction_name2\t$GeneID2\t$geneName2\tNA\t$Sub_interaction_type\tTarget-RBP\tHg19\t$pancancerNum\tNA\tstarBase";
                unless (exists $hash2{$output}){
                    $hash2{$output} =1;
                    print $O1 "$output\n";
                }
            }
            else{
                my $pancancerNum = $f[17];
                my $output = "$chr1\t$start1\t$end1\t$Strand1\t$Type1\t$interaction_name1\t$geneID1\t$geneName1\tNA\t$chr2\t$start2\t$end2\t$Strand2\t$Type2\t$interaction_name2\t$GeneID2\t$geneName2\tNA\t$Sub_interaction_type\tTarget-RBP\tHg19\t$pancancerNum\tNA\tstarBase";
                unless (exists $hash2{$output}){
                    $hash2{$output} =1;
                    print $O1 "$output\n";
                }
            }
        }
    }
}
