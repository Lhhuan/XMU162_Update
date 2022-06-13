# 统计("TSS_TSS","ENH_ENH","TSS_ENH"),enhancer-target 被hotspot cover 的比例，得../output/edges_annotation/12_cover_enhancer_target_type_ratio.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my (%hash1,%hash2,%hash3,%hash4);

my $fo1 = "../output/edges_annotation/12_cover_enhancer_target_type_ratio.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Total_num\tcover_num\tratio\ttype\n";
my @types=("TSS_TSS","ENH_ENH","TSS_ENH");
foreach my $type(@types){

    my $f1= "/home/huanhuan/project/link_database/OncoBase/output/hg38/01_${type}_sorted.bed.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $f2= "../output/edges_annotation/success_${type}_ori.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I1>)
    {
        chomp;
        my @f=split/\t/;
        my $info =$f[-1];
        my @t= split/\;/,$info;
        my $chr1=$t[0];
        my $start1=$t[1];
        my $end1=$t[2];
        my $chr2 =$t[3];
        my $start2 =$t[4];
        my $end2 =$t[5];
        my $gene1 =$t[7];
        my $gene2=$t[8];
        $gene1 =~s/\..*//g;
        $gene2 =~s/\..*//g;
        my $k1 =join("\t",@t[0..5],$gene1,$gene2);
        $hash1{$k1}=1;
    }

    while(<$I2>)
    {
        chomp;
        my @t= split/\;/;
        my $chr1=$t[0];
        my $start1=$t[1];
        my $end1=$t[2];
        my $chr2 =$t[3];
        my $start2 =$t[4];
        my $end2 =$t[5];
        my $gene1 =$t[7];
        my $gene2=$t[8];
        $gene1 =~s/\..*//g;
        $gene2 =~s/\..*//g;
        my $k2 =join("\t",@t[0..5],$gene1,$gene2);
        $hash2{$k2}=1;
    }
    my $num1= keys %hash1;
    my $num2 =keys %hash2;
    my $ratio =$num2/$num1;
    print $O1 "$num1\t$num2\t$ratio\t$type\n";
}

my $f3= "/home/huanhuan/project/link_database/JEME/output/hg38/01_merge_enhancer_target_sample_sorted.bed.gz";
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
my $f4= "../output/nodes_annotation/JEME_enhancer.bed.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

while(<$I3>)
{
    chomp;
    my @f=split/\t/;
    my $k3 =join("\t",@f[0..3]);
    $hash3{$k3}=1;
}

while(<$I4>)
{
    chomp;
    my @f=split/\t/;
    my $k4 =join("\t",@f[4..7]);
    $hash4{$k4}=1;
}

my $num3= keys %hash3;
my $num4 =keys %hash4;
my $ratio =$num4/$num3;
print $O1 "$num3\t$num4\t$ratio\tenhancer_target\n";