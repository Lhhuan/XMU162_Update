#将../raw_data/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt 和../raw_data/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt normalized 成../normalized/01_SEEQTL_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;




my $fo1 = "../normalized/01_SEEQTL_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tP_value\tBeta\tCis_or_Trans\tTissue\tSource";
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
            my $Variant_id=$f[];
            my $Rs_id=$f[];
            my $Variant_Chr=$f[];
            my $Variant_Pos=$f[];
            my $Variant_Ref=$f[];
            my $Variant_Alt=$f[];
            my $Maf=$f[];
            my $Gene_Version=$f[];
            my $Gene_id=$f[];
            my $Gene_name=$f[];
            my $EntrezGeneID=$f[];
            my $Gene_Chr=$f[];
            my $Gene_Start=$f[];
            my $Gene_End=$f[];
            my $Slope=$f[];
            my $Slope_se=$f[];
            my $P_value=$f[];
            my $Beta=$f[];
            my $Cis_or_Trans=$f[];
            my $Tissue=$f[];
            my $Source=$f[];
            my $header = "$Variant_id\t$Rs_id\t$Variant_Chr\t$Variant_Pos\t$Variant_Ref\t$Variant_Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End\t$Slope\t$Slope_se\t$P_value\t$Beta\t$Cis_or_Trans\t$Tissue\t$Source";
        }
    }
}
