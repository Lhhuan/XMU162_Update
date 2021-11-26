#normalize "../raw_data/interactions.tsv.gz", 得../normalized/01_interactions_normalized.txt;
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../raw_data/interactions.tsv.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../normalized/01_interactions_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# my $fo2 = "../output/08_ENCORI_hg19_RNA-RNA_normalized_fail.txt"; #得没有匹配的mirna文件
# open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my (%hash1,%hash2,%hash3);

while(<$I1>) #在06中用end - start确定该文件为0-based
{
    chomp;
    unless(/^ncrna-name/){
        my @f= split/\t/;
        my $ncrna_name1 = $f[0];
        my $ncrna_name2 =$f[3];
        my $Type1 = "NcRNA";
        my $GeneID2 = "NA";
        my $geneName2 ="NA";
        my $GeneID1 = "NA";
        my $geneName1 ="NA";
        my $Type2 = "NcRNA";
        my $chr1 = "NA";
        my $start1 = "NA";
        my $end1 = "NA";
        my $Strand1 = "NA";
        my $interaction_name1 = $ncrna_name1;
        my $chr2 = "NA";
        my $start2 = "NA";
        my $end2 = "NA";
        # my $ddd = $end2 -$start2 ;
        # print "$ddd\n";
        my $Strand2 = "NA";
        my $interaction_name2 = $ncrna_name2;
        my $CellLine_Tissue ="NA";
        my $Sub_interaction_type = "NcRNA-NcRNA";;
        my $output = "$chr1\t$start1\t$end1\t$Strand1\t$Type1\t$interaction_name1\t$GeneID1\t$geneName1\tNA\t$chr2\t$start2\t$end2\t$Strand2\t$Type2\t$interaction_name2\t$GeneID2\t$geneName2\tNA\t$Sub_interaction_type\tRNA-RNA\tNA\t$CellLine_Tissue\tNA\tArena-Idb";
        unless (exists $hash2{$output}){
            $hash2{$output} =1;
            print $O1 "$output\n";
        }
    }
}
