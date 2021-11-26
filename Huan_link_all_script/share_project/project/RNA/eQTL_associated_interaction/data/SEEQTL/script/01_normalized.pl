#将../raw_data/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt 和../raw_data/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt normalized 成../normalized/01_SEEQTL_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;




my $fo1 = "../normalized/01_SEEQTL_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";

my @file_types = ("cis","trans");
foreach my $file_type(@file_types){
    my $f1 = "../raw_data/eQTL_Qvalue_cutoff_hapmap3_${file_type}_hg19.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";

    my(%hash1,%hash2,%hash3,%hash4);
    while(<$I1>){
        chomp;
        unless(/^SNP/){
            my @f =split/\t/;
            my $Variant_id=$f[0];
            my $Rs_id=$f[0];
            my $Variant_Chr=$f[1];
            my $bp_SNP =$f[2];
            my $Variant_Pos=$bp_SNP +1; #此表格中提供的position比dbSNP中报道的小1,所以加1
            my $Variant_Ref="NA";
            my $Variant_Alt="NA";
            my $Effect_allele="NA";
            my $Other_allele="NA";
            my $Maf="NA";
            my $Gene_Version="Hg19";
            my $Gene_id="NA";
            my $Gene_name=$f[-2];
            my $EntrezGeneID=$f[-3];
            my $Gene_Chr="NA";
            my $Gene_Start="NA";
            my $Gene_End="NA";
            my $Slope="NA";
            my $Slope_se="NA";
            my $P_value="NA";
            my $Q_value = $f[-1];
            my $Beta="NA";
            my $Cis_or_Trans=$file_type;
            my $Tissue="NA";
            my $Source="SEEQTL";
            my $Effct_sise ="NA";
            my $Zscore ="NA";
            my $output = "$Variant_id\t$Rs_id\t$Variant_Chr\t$Variant_Pos\t$Variant_Ref\t$Variant_Alt\t$Effect_allele\t$Other_allele\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue\t$Source";
            print $O1 "$output\n";
        }
    }
}
