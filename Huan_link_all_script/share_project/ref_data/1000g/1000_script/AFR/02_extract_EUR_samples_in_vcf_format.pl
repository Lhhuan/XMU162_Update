#将/home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 中提取出EUR 的sample info，得1kg.phase3.v5.shapeit2.eur.hg19.chr${i}.vcf.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager; #多线程并行


my $pm = Parallel::ForkManager->new(22); ## 设置最大的线程数目
my @chr_numbers;
for (my $i=1;$i<23;$i++){ 
    push @chr_numbers,$i;
}

foreach my $i(@chr_numbers){
    my $pid = $pm->start and next; #开始多线程
    my $command = "bcftools view -S EUR_sample_list.txt /home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o 1kg.phase3.v5.shapeit2.eur.hg19.chr${i}.vcf.gz -O z";
    system $command;
    print "$command\n";
    $pm->finish;  #多线程结束
}