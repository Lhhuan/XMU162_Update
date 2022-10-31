 #将"../01_need_cell.txt" 从已有的细胞信息中对应选出来，得02_exists_cell_line_info.txt，没出现的cell得02_NO_info_cell_line.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/04_all_cell_line_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $f3 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_existing_ccle_cell_line_info.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
my $f4 = "../01_need_cell.txt";
open my $I4, '<', $f4 or die "$0 : failed to open input file '$f4' : $!\n"; 

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#-------------------
my $fo1 = "./02_exists_cell_line_info.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

print $O1 "Cell_line\tdisease\n";
my $fo2 = "./02_NO_info_cell_line.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my %hash1;
my %hash2;

while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Cell_line/){
        my $CELL_Line =$f[0];
        my $Disease =$f[1];
        # if ($Disease =~/\bNO\b/){
            $hash1{$CELL_Line}=$Disease;
        # }
    }
}

while(<$I2>){
    chomp;
    my @f = split/\t/;
    unless(/^Cell_line/){
        my $CELL_Line =$f[0];
        my $Disease =$f[1];
        # if ($Disease =~/\bNO\b/){
            $hash1{$CELL_Line}=$Disease;
        # }
    }
}

while(<$I3>){
    chomp;
    my @f = split/\t/;
    unless(/^Cell_line/){
        my $CELL_Line =$f[0];
        my $Disease ="Cancer";
        # if ($Disease =~/\bNO\b/){
            $hash1{$CELL_Line}=$Disease;
        # }
    }
}

while(<$I4>){
    chomp;
    my $cell =$_;
    if(exists $hash1{$cell}){
        my $v = $hash1{$cell};
        print $O1 "$cell\t$v\n";
    }
    else{
        print $O2 "$cell\n";
    }
}


# my %hash4;
# for (my $i=1;$i<23;$i++){
#     my $k = "chr${i}";
#     $hash4{$k}=1;
# }

# while(<$I2>){
#     chomp;
#     my @f = split/\s+/;
#     my $link = $f[2];
#     my @t = split/\//,$link;
#     # $link =~ s///g;
#     # print  "$t[7]\n"; 
#     my $org_file = $t[7];
#     my $file_name = $org_file;
#     $file_name =~ s/wgEncodeUwTfbs//g;   
#     $file_name =~ s/Ctcf.*//g;   
#     # wgEncodeUwTfbsA549CtcfStdPkRep2.narrowPeak.gz
#     if(exists $hash1{$file_name }){
#         my $f3 = "../data/$org_file";
#         print  "$f3\n";
#         # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
#         open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
#         while(<$I3>)
#         {
#             chomp;
#             # unless(/^chr/){
#             my @f = split/\s+/;
#             my $chr = $f[0];
#             my $start =$f[1];
#             my $end =$f[2];
#             my $output1 = "$chr\t$start\t$end";
#             if (exists $hash4{$chr}){
#                 unless(exists $hash2{$output1}){
#                     $hash2{$output1}=1;
#                     print $O1 "$output1\n";
#                 }
#             }
#         }
#     }
    
    
# }