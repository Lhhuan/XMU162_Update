#"../output/07_all_predict_top0.1_need_left_overlap_hi_c.bed.gz" header 为$h_chr\t$h_start\t$h_end\t$egene\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end\thi-c_value
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

my $f1 = "../output/07_groundtruth_chr1_left_overlap_hi_c.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# my $f2 = "../output/05_all_predict_top0.1.bed.gz";
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/07_groundtruth_chr1_gene_hotspot.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/07_hic_overlap_groundtruth_chr1_left_right_left.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my (%hash1,%hash2);


while(<$I1>)
{
    chomp;
    unless(/^hotspot_id/){
        my @f=split/\t/;
        my $h_chr =$f[0];
        my $h_start =$f[1];
        my $h_end=$f[2];
        my $egene =$f[3];
        my $gene_chr =$f[4];
        $gene_chr = "chr${gene_chr}";
        my $gene_start =$f[5];
        my $gene_end =$f[6];
        my $segment1_chr =$f[7];
        my $segment1_start =$f[8];
        my $segment1_end =$f[9];
        my $segment2_chr =$f[10];
        my $segment2_start =$f[11];
        my $segment2_end = $f[12];
        my $output1 = "$gene_chr\t$gene_start\t$gene_end\t$egene\t$h_chr\t$h_start\t$h_end";
        my $output2= "$segment2_chr\t$segment2_start\t$segment2_end\t$segment1_chr\t$segment1_start\t$segment1_end";
        unless(exists $hash1{$output1}){
            $hash1{$output1}=1;
            print $O1 "$output1\n";
            # print  "$output1\n";
        }
        unless(exists $hash2{$output2}){
            $hash2{$output2}=1;
            print $O2 "$output2\n";
        }
    }
}

close($O1);
close($O2);

my $sort_fo1 = "../output/07_groundtruth_chr1_gene_hotspot_sorted.bed.gz";
my $sort_fo2 = "../output/07_hic_overlap_groundtruth_chr1_left_right_left_sorted.bed.gz";
my $rightoverlp = "../output/07_groundtruth_chr1_right_overlap_hi_c.bed.gz";
system "zless $fo1 |sort -k1,1 -k2,2n |gzip >$sort_fo1";
system "zless $fo2 |sort -k1,1 -k2,2n |gzip >$sort_fo2";
system "bedtools intersect -a $sort_fo1 -b $sort_fo2 -wo |gzip > $rightoverlp";





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