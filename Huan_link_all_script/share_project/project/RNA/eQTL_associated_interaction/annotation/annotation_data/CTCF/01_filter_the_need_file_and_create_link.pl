#将files.txt中ctcf对于的narrow peak选出来，对于同一样本有重复的，筛选出文件大的哪个，对于并生成download link,得download.sh,和need_file.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "files.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "need_file.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "download.sh";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

my (%hash1,%hash2,%hash3,%hash4);

while(<$I1>)
{
    chomp;
    my @f = split/\s+/;
    my $file_name = $f[0];
    my $size =$f[-1];
    $size =~ s/size=//g;
    if ($file_name =~/Ctcf/ && $file_name=~/narrowPeak/ ){
    # if ($file_name =~/Ctcf/ && $file_name=~/narrowPeak/ && $_=~/origAssembly=hg19/ ){
        # print "$file_name\n";
        my $file = $file_name;
        $file =~ s/Rep.*+$//g; #从Rep 到最后替换成空
        my $k1 = "$file\t$size";
        # print "$k1\n";
        # push @{$hash1{$k1}},$_;
        $hash1{$k1}=$_;#为hash3下面取整行数据，做准备
        push @{$hash2{$file}},$size;
    }
}


foreach my $file(sort keys %hash2){ #选出每个sample 对应的最大的Rep
    my @sizes = @{$hash2{$file}};
    my $number =@sizes;
    if ($number >1){
        my $all_size = join("\t",@sizes);
        print "$all_size\n";
        if ($all_size =~/M/){ #size 是1MB的那几个先挑选出来(每个数组内只要一个MB,其他是kb)
            # print "$all_size\n";
            foreach my $size(@sizes){
                if ($size =~/M/){
                    my $k3 = "$file\t$size";
                    $hash3{$k3}=1                   
                }
            }
        }
        else{
            my @sizes = sort {$b cmp $a} @sizes; #对数组内元素按照字符串降序排序
            my $size = $sizes[0];#取最大值
            my $k3 = "$file\t$size";
            $hash3{$k3}=1
        }
    }
    else{ #只有一个size
        my $size = $sizes[0];
        my $k3 = "$file\t$size";
        $hash3{$k3}=1
    }
}


foreach my $k3 (sort keys %hash3){
    if (exists $hash1{$k3}){
        my $v = $hash1{$k3};
        print $O1 "$v\n";
        my @f = split/\s+/,$v;
        my $file_name = $f[0];
        my $command = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/${file_name}"; 
        print $O2 "wget -c $command\n";
    }
}
