 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../output/Whole_Blood_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
    my $dir = "../../data/GTEx_Analysis_v8_eQTL_hg19";
    my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
    my $tissue = "Whole_Blood";
    my $f1 ="${dir}/${tissue}${suffix}";
    # my $f1 = "/share/data0/QTLbase/2020_12_2_QC/${QTL}.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $f2 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    #---------------------------
    # my $output_dir = "../output/ALL_${QTL}";
    # if(-e $output_dir){
    #     print "${output_dir}\texist\n";
    # }
    # else{
    #     system "mkdir -p $output_dir";
    # }
    #---------------
    my $fo1 = "../../output/Whole_Blood_cis_eQTL_1kg_Completion.txt.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "SNP_chr\tSNP_pos\tPvalue\n";

    my %hash1;
    while(<$I1>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $Pvalue =$f[-6];
            my $k = "$SNP_chr\t$SNP_pos";
            $hash1{$k}=1;
            print $O1 "$k\t$Pvalue\n";
        }
    }


    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
            }
        }
    }
    # $pm->finish;  #多线程结束
# }

