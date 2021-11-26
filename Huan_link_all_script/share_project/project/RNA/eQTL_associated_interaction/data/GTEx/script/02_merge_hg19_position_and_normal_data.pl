#将../output/01_unique_variant-id_position_hg19.vcf和../output/01_normal_format.txt.gz merge 在一起，得../normalized/02_GTEx_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/01_unique_variant-id_position_hg19.vcf";
my $f2 = "../output/01_normal_format.txt.gz";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../normalized/02_GTEx_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";

my(%hash1,%hash2,%hash3,%hash4);
while(<$I1>){
    chomp;
    unless(/^#/){
        my @f =split/\s+/;
        my $chr= $f[0];
        my $pos = $f[1];
        my $ID = $f[2];
        $hash1{$ID}= $pos;
    }
}

while(<$I2>){
    chomp;
    my @f= split/\t/;
    unless(/^Variant_id/){
        my @f = split/\t/;
        my $Variant_id =$f[0];
        my $Rs_id =$f[1];
        my $Chr =$f[2];
        my $Pos =$f[3];
        my $ref =$f[4];
        my $alt = $f[5];
        my $Effect_allele=$alt;#网址记录：the eQTL effect allele is the ALT allele
        my $Other_allele="NA";
        my $Maf =$f[6];
        my $Gene_Version =$f[7];
        my $Gene_id =$f[8];
        my $Gene_name = $f[9];
        my $EntrezGeneID ="NA";
        my $Gene_Chr ="NA";
        my $Gene_Start ="NA";
        my $Gene_End ="NA";
        my $gene_info ="$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End";
        my $other_info = join("\t",@f[10..$#f]); #第11个元素到最后
        if (exists $hash1{$Variant_id}){ #hg38转hg19成功的id
            my $pos_hg19 = $hash1{$Variant_id};
            $Gene_Version = "Hg19";
            my $Variant_id_hg19 = "${Chr}_${pos_hg19}_${ref}_${alt}_b19";
            my $output = "$Variant_id_hg19\t$Rs_id\t$Chr\t$pos_hg19\t$ref\t$alt\t$Effect_allele\t$Other_allele\t$Maf\t$Gene_Version\t$gene_info\t$other_info\tGTEx";
            print $O1 "$output\n";
        }
        else{#hg38转hg19失败的id
            my $varint_info = join ("\t",@f[0..5],"NA","NA",@f[6,7]);
            print $O1 "$varint_info\t$gene_info\t$other_info\tGTEx\n";
        }
    }
}

