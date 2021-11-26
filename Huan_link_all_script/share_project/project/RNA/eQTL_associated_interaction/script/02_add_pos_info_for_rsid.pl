#将用"~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz" 为 ../output/01_all_eQTL_associated_interaction.txt.gz填充position和ref 和alt的信息
#得../output/02_add_pos_info_for_rsid.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/01_all_eQTL_associated_interaction.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/02_add_pos_info_for_rsid.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";


while(<$I1>)
{
    chomp;
    unless(/^Chr/){
        my @f =split/\t/;
        my $Chr = $f[0];
        my $Pos =$f[1];
        my $RS_id =$f[2];
        my $Ref =$f[3];
        my $Alt =$f[4];
        my $v = "$Chr\t$Pos\t$Ref\t$Alt";
        $hash1{$RS_id}=$v;
    }
}


while(<$I2>)
{
    chomp;
    unless(/^Variant_id/){
        my @f =split/\t/;
        my $Variant_id = $f[0];
        my $Rs_id =$f[1];
        my $Variant_Chr =$f[2];
        my $Variant_Pos =$f[3];
        my $Variant_Ref =$f[4];
        my $Variant_Alt =$f[5];
        my $other_info = join("\t",@f[6..$#f]);
        my $variant_info = join("\t",@f[2..5]);
        if ($variant_info =~/NA/){
            if ($hash1{$Variant_id}){
                my $v = $hash1{$Variant_id};
                my @t = split/\t/,$v;
                my $v_chr =$t[0];
                my $v_pos =$t[1];
                my $v_ref =$t[2];
                my $v_alt =$t[3];
                if ($Variant_Alt =~/NA/){ #如果Variant_Alt 为NA,则原表格中的chr,pos,ref,alt都被替换
                    my $output = "$Variant_id\t$Rs_id\t$v\t$other_info";
                    print $O1 "$output\n";
                }
                else{#如果Variant_Alt 不为NA，则保留原来的alt,因为同一位置不同alt的snp对应的是同一个rsid,此处只替换$Chr\t$Pos\t$Ref
                    my $output= "$Variant_id\t$Rs_id\t$v_chr\t$v_pos\t$v_ref\t$Variant_Alt\t$other_info";
                    print $O1 "$output\n";
                }
            }
            else{ #chr,pos,ref,alt中NA,但是不在key 不在hash1。原文输出
               print $O1 "$_\n";
            }
        }
        else{ # chr,pos,ref,alt均不为NA
            print $O1 "$_\n";
        }
    }
}

close $O1;

    