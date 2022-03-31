#"../output/07_random_top0.1_need_left_overlap_hi_c.bed.gz" header 为$h_chr\t$h_start\t$h_end\t$egene\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end\thi-c_value
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/script/groundtruth/output/hi-c_top25/01_Hi_C_result_chr1_top25.bed.gz";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "../output/07_random_chr1.bed.gz";
my $sort_f2 = "../output/07_random_top0.1_need_sorted.bed.gz";
my $leftoverlp = "../output/07_random_top0.1_need_left_overlap_hi_c.bed.gz";
system "zless $fo1 |sort -k1,1 -k2,2n |gzip >$sort_f2";

system "bedtools intersect -a $sort_f2 -b $f1 -wo |gzip > $leftoverlp";





# my $fo1 = "../output/05_trans_predict.txt.gz";
# open my $O1, "| gzip >$fo1" or die $!;
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# my $fo2 = "../output/05_cis_predict.txt.gz";
# open my $O2, "| gzip >$fo2" or die $!;

# my $header ="hotspot_id\tgene_id\tsmiliarity\th_chr\th_start\th_end\tgene\tgene_chr\tgene_start\tgene_end";
# print $O1 "$header\n";
# print $O2 "$header\n";

# my (%hash1,%hash2,%hash3,%hash4);

# while(<$I1>)
# {
#     chomp;
#     unless(/^segment1/){
#         my @f=split/\t/;
#         my $segment1 =$f[0];
#         my $segment2 =$f[1];
#         my $value=$f[2];
#         my @t1=split/_/,$segment1;
#         my $chr1 = $t1[0];
#         my $start1 = $t1[1];
#         my $end1 = $t1[2];
#         my @t2=split/_/,$segment2;
#         my $chr2 = $t2[0];
#         my $start2 = $t2[1];
#         my $end2 = $t2[2];
#         my $output = 
#     }
# }