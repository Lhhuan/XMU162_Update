#将/home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 中提取出SAS 的sample info，得1kg.phase3.v5.shapeit2.sas.hg19.chr${i}.vcf.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
# use Parallel::ForkManager; #多线程并行


# my $pm = Parallel::ForkManager->new(22); ## 设置最大的线程数目
my @chr_number=(1..22);
my @chr_numbers = sort { $b <=> $a } @chr_number;


# foreach my $chr(@chr_numbers){
#     print "$chr\n";
# }



foreach my $i(@chr_numbers){
    my $command = "bcftools view -S EUR_sample_list.txt ../CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz -o EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz";
    print "$command\n";
    system $command;
    print "$i\n";
    # $pm->finish;  #多线程结束
}