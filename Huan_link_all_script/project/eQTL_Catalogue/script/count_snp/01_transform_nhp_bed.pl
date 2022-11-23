#将../../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz 转化为bed格式，得./output/01_nhp.bed.gz,排序得./output/01_nhp_sorted.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

# my $pm = Parallel::ForkManager->new(30); ## 设置最大的线程数目

my $fo2 = "./output/01_nhp.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# print $O2 "SNP_chr\tSNP_pos\tPvalue\n";

my $f1 = "../../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    unless(/^emplambda/){
        my $chr =$f[2];
        my $CHR = "chr${chr}";
        my $POS = $f[1];
        my $start =$POS-1; #避免文件出现科学计数法(e+)
        my $end = $POS+1-1;#避免文件出现科学计数法(e+)
        print $O2 "$CHR\t$start\t$end\n";
    }
}

close($I1);
close($O2);
system "zless $fo2 |sort -k1,1 -k2,2n |gzip >./output/01_nhp_sorted.bed.gz"