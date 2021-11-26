#将用"~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz" 为 ../output/02_add_pos_info_for_rsid.txt.gz填充position和ref 和alt的信息 #得../output/03_add_rsid.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "~/ref_data/dbSNP/01_extract_b152_vcf.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/01_all_eQTL_associated_interaction.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/03_add_rsid.txt.gz";
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
        my @alts = split/\,/,$Alt;
        foreach my $alt(@alts){
            my $key = "$Chr\t$Pos\t$Ref\t$alt";
            $hash1{$key}=$RS_id;
        }
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
        my $Variant = join ("\t",@f[2..5]);
        my $other_info = join("\t",@f[6..$#f]);
        if ($Rs_id =~/rs/){ #有rs id 的就不再寻找
            print $O1 "$_\n";
        }
        elsif($Rs_id =~/NA/){ #rsID =NA
            my $Variant_info = $Variant;
            $Variant_info=~s/chr//g;
            if (exists $hash1{$Variant_info}){ #可以在hash1中找到rs id
                my $id = $hash1{$Variant_info};
                my $output = "$Variant_id\t$id\t$Variant\t$other_info";
                print $O1 "$output\n";
            }
            else{ #不能在hash1中找到rs id,原文输出
                print $O1 "$_\n";
            }
        }
        else{ #rsID =1111； 没有rsid,需要把rsid 变为rs1111
            $Rs_id ="rs${Rs_id}";
            my $output = "$Variant_id\t$Rs_id\t$Variant\t$other_info";
            print $O1 "$output\n";
        }
    }
}

close $O1;

    