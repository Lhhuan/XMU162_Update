#将metadata.tsv中选一个数据最大的出来，对于同一样本有重复的，筛选出文件大的哪个，对于并生成download link,得download.sh,和need_metadata.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "metadata.tsv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "need_metadata.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "download.sh";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
print $O1 "Biosample_term_name\tExperiment_target\tLab\tsize\tFile_accession\tLab\tFile_download_URL\n";
my (%hash1,%hash2,%hash3,%hash4);

while(<$I1>)
{
    chomp;
    unless(/^File/){
        my @f = split/\t/;
        my $File_accession = $f[0];
        my $Output_type =$f[4];
        my $File_assembly = $f[5];
        my $Biosample_term_name = $f[9];
        my $Experiment_target =$f[21];
        my $size = $f[42];
        my $Lab = $f[43];
        my $File_download_URL = $f[46];
        $File_download_URL =~s/\@//g;
        if($Output_type =~/^\bpeaks\b$/ && $File_assembly =~/hg19/){
            # print "$File_download_URL\n";
            my $k1 = "$Biosample_term_name\t$Experiment_target\t$Lab\t$size";
            my $v1 = "$File_accession\t$Lab\t$File_download_URL";
            $hash1{$k1}=$v1;
            my $k2 = "$Biosample_term_name\t$Experiment_target\t$Lab";
            push @{$hash2{$k2}},$size;
        }
    }
}

foreach my  $k2 (sort keys %hash2){
    my @sizes = @{$hash2{$k2}};
    my $number = @sizes;
    if ($number >1){
        my @sorted_sizes = sort {$b <=> $a} @sizes; #对数组内元素按照数字降序排序
        my $size= $sorted_sizes[0]; #一组中选size最大的哪个file
        my $k3 ="$k2\t$size";
        $hash3{$k3}=1;
    }
}

foreach my $k3(sort keys %hash3){
    if (exists $hash1{$k3}){
        my $v = $hash1{$k3};
        print $O1 "$k3\t$v\n";
        my @t =split/\t/,$v;
        my $link = $t[-1];
        print $O2 "wget -c $link\n";
    }
}



# foreach my $file(sort keys %hash2){ #选出每个sample 对应的最大的Rep
#     my @sizes = @{$hash2{$file}};
#     my $number =@sizes;
#     if ($number >1){
#         my $all_size = join("\t",@sizes);
#         print "$all_size\n";
#         if ($all_size =~/M/){ #size 是1MB的那几个先挑选出来(每个数组内只要一个MB,其他是kb)
#             # print "$all_size\n";
#             foreach my $size(@sizes){
#                 if ($size =~/M/){
#                     my $k3 = "$file\t$size";
#                     $hash3{$k3}=1                   
#                 }
#             }
#         }
#         else{
#             my @sizes = sort {$b cmp $a} @sizes; #对数组内元素按照字符串降序排序
#             my $size = $sizes[0];#取最大值
#             my $k3 = "$file\t$size";
#             $hash3{$k3}=1
#         }
#     }
#     else{ #只有一个size
#         my $size = $sizes[0];
#         my $k3 = "$file\t$size";
#         $hash3{$k3}=1
#     }
# }


# foreach my $k3 (sort keys %hash3){
#     if (exists $hash1{$k3}){
#         my $v = $hash1{$k3};
#         print $O1 "$v\n";
#         my @f = split/\s+/,$v;
#         my $file_name = $f[0];
#         my $command = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/${file_name}"; 
#         print $O2 "wget -c $command\n";
#     }
# }
