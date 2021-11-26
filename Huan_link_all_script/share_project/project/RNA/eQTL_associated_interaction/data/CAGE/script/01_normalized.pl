#将../raw_data/cage_association_summary_statistics.txt normalized 成../normalized/01_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;




my $fo1 = "../normalized/01_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";


my $f1 = "../raw_data/cage_association_summary_statistics.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";

my(%hash1,%hash2,%hash3,%hash4);
while(<$I1>){
    chomp;
    unless(/^SNP/){
        my @f =split/\t/;
        my $Variant_id=$f[0];
        my $Rs_id=$f[0];
        my $Variant_Chr=$f[1];
        my $Variant_Pos=$f[2];
        my $Variant_Ref="NA";
        my $Variant_Alt="NA";
        my $Effect_allele=$f[3];
        my $Other_allele=$f[4];
        my $Maf="NA";
        my $Gene_Version="Hg19";
        my $Gene_id=$f[-1];
        my $Gene_name=$f[-2];
        my $EntrezGeneID="NA";
        my $Gene_Chr="NA";
        my $pos ="NA";
        my $Gene_Start = "NA";
        my $Gene_End= "NA";
        my $Slope="NA";
        my $Slope_se= $f[8];
        my $P_value= $f[10];
        my $Q_value = "NA";
        my $Beta= $f[7];
        my $Cis_or_Trans= "NA";
        my $Tissue="Peripheral Blood";
        my $Source="CAGE";
        my $Effct_sise ="NA";
        my $Zscore ="NA";
        my $output = "$Variant_id\t$Rs_id\t$Variant_Chr\t$Variant_Pos\t$Variant_Ref\t$Variant_Alt\t$Effect_allele\t$Other_allele\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue\t$Source";
        print $O1 "$output\n";
    }
}
