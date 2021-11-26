# 统计"$output_dir/${tissue}_random_select_result_${i}.txt.gz" 中，hit住eQTL的数目，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/count_random_hit_eQTL_ratio.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/count_random_hit_eQTL_ratio.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Tissue\tRandom_number\torg_eqtl_count\trandom_hit_eQTL_count\tRandom_hit_eQTL_ratio\n";


my $dir = "../../data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
my @tissues = ("Whole_Blood", "Lung");

my %hash1;
for(my $i=1; $i<23;$i++){
    $hash1{$i}=1;
}


foreach my $tissue (@tissues){
    my $f1 ="${dir}/${tissue}${suffix}";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my(%hash2,%hash3,%hash4);
    while(<$I1>){
        chomp;
        my @f = split/\t/;
        unless(/^variant_id/){
            my $variant_id = $f[0];
            my $chr = $f[1];
            my $pos = $f[2];
            my $Pvalue =$f[-6];
            if (exists $hash1{$chr}){
                if($Pvalue <5e-8){ #----significant qtl
                    my $k = "$chr\t$pos";
                    # print "$_\n";
                    $hash2{$k}=1;
                    print "$chr\n";
                }
            }
        }
    }
    my $org_eqtl_count = keys %hash2;
    print "$tissue\t$org_eqtl_count\n";
}