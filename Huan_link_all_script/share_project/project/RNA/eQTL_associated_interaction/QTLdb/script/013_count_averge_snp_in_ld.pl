# 计算../output/012_random_sample_snp_for_ld_block_r2_${i}.tags.list $i =0.5,0.7,0.8 中ld中平均包含的snp的个数得../output/013_count_averge_snp_in_ld.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "../output/013_count_averge_snp_in_ld.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
print $O1 "r2\taverage_number_of_snp\n";

my @ccc= (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8);
    foreach my $i(@ccc){
    my $f1 = "../output/012_random_sample_snp_for_ld_block_r2_${i}.tags.list";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my %hash1;
    my @snps_count;

    while(<$I1>)
    {
        chomp;
        unless(/^SNP/){
            my @f = split/\s+/;
            my $SNP =$f[1];
            my $TAGS =$f[-1]; 
            # print "$SNP\n";
            # print "$SNP\t$TAGS\n";
            if($TAGS =~/NONE/){ #   1 snp in the ld
            #    print "$SNP\t$TAGS\n"; 
            push @snps_count,1;
            }
            else{
                my @t= split/\|/,$TAGS;
                push @t,$SNP;
                my %hash2;
                my @uniq_t = grep { ++$hash2{ $_ } < 2; } @t;
                my $number = @uniq_t;
                push @snps_count,$number;
            }
        }
    }
    my @sorted_snps_count = sort {$b <=> $a} @snps_count;
    my @final = @sorted_snps_count[1..$#sorted_snps_count]; #去掉了最大值
    my $sum= sum(@final);
    my $length = @final;
    my $average = $sum/$length;
    print $O1 "$i\t$average\n";
    # print "@final\n";
    # print "$length\n";
}