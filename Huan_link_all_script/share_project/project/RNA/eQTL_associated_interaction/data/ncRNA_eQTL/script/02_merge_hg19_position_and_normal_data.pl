#将../output/01_unique_gene_position_hg19.bed和../output/01_normal_format.txt.gz merge 在一起，得../normalized/02_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/01_unique_gene_position_hg19.bed";
my $f2 = "../output/01_normal_format.txt.gz";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../normalized/02_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";


my(%hash1,%hash2,%hash3,%hash4);
while(<$I1>){
    chomp;
    unless(/^#/){
        my @f =split/\s+/;
        my $chr= $f[0];
        my $start = $f[1];
        my $end= $f[2];
        my $gene_id =$f[3];
        my $v= "$chr\t$start\t$end";
        $hash1{$gene_id}= $v;
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
        my $Effect_allele=$ref;
        my $Other_allele=$alt;
        my $Maf =$f[6];
        my $Gene_Version ="Hg19";
        my $Gene_id =$f[8];
        my $Gene_name = $f[9];
        my $EntrezGeneID ="NA";
        # my $Gene_Chr ="NA";
        # my $Gene_Start ="NA";
        # my $Gene_End ="NA";
        my $other_info = join("\t",@f[10..$#f]); #第11个元素到最后
        if (exists $hash1{$Gene_id}){ #hg38转hg19成功的gene id
            my $gene_pos_hg19 = $hash1{$Gene_id};
            my @t = split/\t/,$gene_pos_hg19;
            my $Gene_Chr =$t[0];
            my $Gene_Start = $t[1];
            my $Gene_End =$t[2];
            my $gene_info ="$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End";
            my $output = "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$ref\t$alt\t$Effect_allele\t$Other_allele\t$Maf\t$Gene_Version\t$gene_info\t$other_info\tncRNA_eQTL";
            unless(exists $hash2{$output}){
                $hash2{$output} =1;
                print $O1 "$output\n";
            }
        }
        else{#hg38转hg19失败的id
            my $Gene_Chr ="NA";
            my $Gene_Start ="NA";
            my $Gene_End ="NA";
            my $gene_info ="$Gene_id\t$Gene_name\t$EntrezGeneID\t$Gene_Chr\t$Gene_Start\t$Gene_End";
            my $output = "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$ref\t$alt\t$Effect_allele\t$Other_allele\t$Maf\t$Gene_Version\t$gene_info\t$other_info\tncRNA_eQTL";
            print $O1 "$output\n";
        }
    }
}

