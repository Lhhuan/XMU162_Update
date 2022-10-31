#将files.txt中未经其他处理的bed选出来，对于同一样本有重复的，筛选出文件大的哪个，对于并生成download link,得download.sh,和01_need_file.txt 和对应的细胞文件01_need_cell.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "files.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "01_need_file.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "download.sh";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 = "01_need_cell.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my (%hash1,%hash2,%hash3,%hash4);

while(<$I1>)
{
    chomp;
    my @f = split/\s+/;
    my $file_name = $f[0];
    my $cell =$f[7];
    my $treatment=$f[8];
    my $file_type= $f[-3];
    my $size =$f[-1];
    # if()
    $size =~ s/size=//g;

    # print "$file_name\t$cell\t$treatment\t$size\t$file_type\n";
    # # print "$size\n";
    if ($treatment =~/treatment=None/ && $file_type=~/type=bed/ ){
        my $k1 = "$cell\t$size";
        # print "$size\n";
        $hash1{$k1}=$_;#为hash3下面取整行数据，做准备
        push @{$hash2{$cell}},$size;
    }
}


foreach my $file(sort keys %hash2){ #选出每个sample 对应的最大的Rep
    my @sizes = @{$hash2{$file}};
    my $number =@sizes;
    if ($number >1){
        my @sizes = sort {$b cmp $a} @sizes; #对数组内元素按照字符串降序排序
        my $size = $sizes[0];#取最大值
        my $k3 = "$file\t$size";
        $hash3{$k3}=1
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
        my $command = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/${file_name}"; 
        print $O2 "wget -c $command\n";
        my @t =split/\s+/,$k3;
        my $cell =$t[0];
        $cell =~ s/cell=|;//g;
        print $O3 "$cell\n";
    }
}
