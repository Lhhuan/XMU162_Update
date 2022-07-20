#统计"./output/hotspot_snp.bed.gz" 包含的snps数目，得"./output/02_count_snp_in_hotspot.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

# my $pm = Parallel::ForkManager->new(30); ## 设置最大的线程数目

my $fo2 = "./output/02_count_snp_in_hotspot.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "chr\tstart\tend\tSNP_number\tcenter_snp\n";

my $f1 = "./output/hotspot_snp.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $k = join("\t",@f[0..2]);
    # my $v = join("\t",@f[3..5]);
    my $v= $f[4];
    push @{$hash1{$k}},$v;
}
foreach my $k(keys %hash1){
    my @vs =@{$hash1{$k}};
    my %ha;
    my @vs1=grep{++$ha{$_}<2}@vs;
    my $num =@vs1;
    # print $O2 "$k\t$num\n";
    my @sorted_vs = sort {$a <=> $b} @vs1;
    my $int_n = int($num/2);
    my $remainder_n =$num % 2;
    if($remainder_n>0){ #奇数
        my $center=$int_n+$remainder_n;
        my $v_center = $sorted_vs[$center-1]; #0-based
        print $O2 "$k\t$num\t$v_center\n";
    }
    else{#偶数
        # print "$int_n\t$num\n";
        my $center=$int_n; #下一步，前9后8
        my $v_center = $sorted_vs[$center];
        print $O2 "$k\t$num\t$v_center\n";
    }
}