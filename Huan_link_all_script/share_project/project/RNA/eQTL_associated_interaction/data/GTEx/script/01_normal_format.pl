#将../raw_data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz进行初步normal, 得gene name 和部分RSid从对应的v8.egenes.txt.gz中获得，得文件../output/01_normal_format.txt.gz
#并得unique的varint POS的vcf文件../output/01_unique_variant-id_position_hg38.vcf.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


#-------------------------------------------------------------#获取组织名称
my $dir = "../raw_data/GTEx_Analysis_v8_eQTL";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
foreach my $file(@files){
    my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
    if ($file =~/$suffix/){ 
        $file =~ s/$suffix//g; #提取前缀disease
        my $tissue =$file;
        push @tissues,$tissue;
    }
}
#---------------------------------------------------#对文件进行处理
my $fo1 = "../output/01_normal_format.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/01_unique_variant-id_position_hg38.vcf.gz";
open my $O2, "| gzip >$fo2" or die $!;
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $header = "Variant_id\tRs_id\tChr\tPos\tRef\tAlt\tMaf\tGene_Version\tGene_id\tGene_name\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue";
print $O1 "$header\n";
print $O2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

my(%hash1,%hash2,%hash3,%hash4);
foreach my $tissue(@tissues){
    my $f1 = "${dir}/${tissue}.v8.egenes.txt.gz";
    my $f2 ="${dir}/${tissue}${suffix}";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

    while(<$I1>){
        chomp;
        my @f = split/\t/;
        unless(/^gene_id/){
            my $Gene_id = $f[0];
            my $Gene_name =$f[1];
            my $Variant_id =$f[11];
            my $Rs_id= $f[18];
            my $chr= $f[13];
            my $pos =$f[14];
            my $ref = $f[15];
            my $alt =$f[16];
            $hash1{$Gene_id}=$Gene_name;
            # my $varint_v = "$Rs_id\t$chr\t$pos\t$ref\t$alt";
            $hash2{$Variant_id}=$Rs_id;
            # print "$Variant_id\t$Rs_id\n";
        }
    }
    # print "$f1\n";
    while(<$I2>){
        chomp;
        my @f =split/\t/;
        unless(/^variant_id/){
            my $Variant_id =$f[0];
            my @t = split/\_/,$Variant_id;
            my $Chr =$t[0];
            my $Pos =$t[1];
            my $Ref=$t[2];
            my $Alt =$t[3];
            my $Gene_id = $f[1];
            my $Maf =$f[5];
            my $pval_nominal =$f[6];
            my $P_value =$pval_nominal; 
            my $Slope =$f[7]; #斜率
            my $Slope_se = $f[8];
            my $pval_beta =$f[11];
            my $Beta = $pval_beta;
            my $Gene_Version ="Hg38";
            my $Cis_or_Trans = "cis";
            my $Q_value = "NA";
            my $Effct_sise =$Slope;#The normalized effect size (NES) of the eQTLs is defined as the slope of the linear regression,
            my $Zscore= "NA";
            if (exists $hash1{$Gene_id}){
                my $Gene_name = $hash1{$Gene_id};
                if (exists $hash2{$Variant_id}){
                    my $Rs_id = $hash2{$Variant_id};
                    my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$tissue";
                    print $O1 "$output\n";
                }
                else{
                    my $Rs_id  = "NA"; #这里不再单独去找$Rs_id，最后汇总一起找
                    my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$tissue";
                    print $O1 "$output\n";
                }
            }
            else{
                print "$Gene_id\n";
            }
            #-----------------------构成vcf文件，用于转hg19;
            my $variant_info = "$Chr\t$Pos\t$Variant_id\t$Ref\t$Alt\t.\t.\t.";
            unless (exists $hash3{$variant_info}){
                $hash3{$variant_info} =1;
                print $O2 "$variant_info\n";
            }
        }
    }
    close $I1;
    close $I2;
}
close $O1;
