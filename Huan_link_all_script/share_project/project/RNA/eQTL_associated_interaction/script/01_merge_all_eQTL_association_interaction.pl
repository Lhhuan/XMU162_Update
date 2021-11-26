#将../data/* 下面normalizeda好的文件合在一起,得../output/01_all_eQTL_associated_interaction.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @datasets = ("Blood_eQTL_browser","CAGE","eQTL_Catalog","eQTL_databases","eQTLGen","GTEx","PancanQTL","SEEQTL","ncRNA_eQTL");
foreach my $dataset(@datasets){
    my $path = "../data/${dataset}/normalized";
    my $command = "zcat ${path}/* >> ../output/01_all_eQTL_associated_interaction.txt";
    system "$command";
}

my $f1 = "../output/01_all_eQTL_associated_interaction.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/01_all_eQTL_associated_interaction.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";

while(<$I1>) #去掉重复的header
{
    chomp;
    unless(/^Variant_id/){
        print $O1 "$_\n";

    }
}

close $O1;
system "rm ../output/01_all_eQTL_associated_interaction.txt"; #删掉中间文件，节省空间
    