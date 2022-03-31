#"../output/07_all_predict_top0.1_need_left_overlap_hi_c.bed.gz" header 为$h_chr\t$h_start\t$h_end\t$egene\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end\thi-c_value
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

my $f1 = "../output/07_groundtruth_chr1_right_overlap_hi_c.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../output/07_groundtruth_chr1_left_overlap_hi_c.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/08_groundtruth_chr1_overlap_hi_c.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# my $fo2 = "../output/07_hic_overlap_groundtruth_chr1_left_right_left.bed.gz";
# open my $O2, "| gzip >$fo2" or die $!;
my (%hash1,%hash2);
print $O1 "h_chr\th_start\th_end\tgene_chr\tgene_start\tgene_end\tsegment1_chr\tsegment1_start\tsegment1_end\tsegment2_chr\tsegment2_start\tsegment2_end\tleft_overlap_bp\tright_overlap_bp\n";

while(<$I1>)
{
    chomp;
    unless(/^hotspot_id/){
        my @f=split/\t/;
        my $gene_chr =$f[0];
        my $gene_start =$f[1];
        my $gene_end =$f[2];
        my $egene =$f[3];
        my $h_chr =$f[4];
        my $h_start =$f[5];
        my $h_end=$f[6];
        my $segment2_chr =$f[7];
        my $segment2_start =$f[8];
        my $segment2_end = $f[9];
        my $segment1_chr =$f[10];
        my $segment1_start =$f[11];
        my $segment1_end =$f[12];
        my $overlap_bp = $f[-1];
        my $k= "$h_chr\t$h_start\t$h_end\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end";
        $hash1{$k}=$overlap_bp;
    }
}

#---------------
while(<$I2>)
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
        my $left_overlap_bp = $f[-1];
        # my $output1 = "$gene_chr\t$gene_start\t$gene_end\t$egene\t$h_chr\t$h_start\t$h_end";
        # my $output2= "$segment2_chr\t$segment2_start\t$segment2_end\t$segment1_chr\t$segment1_start\t$segment1_end";
        my $k= "$h_chr\t$h_start\t$h_end\t$gene_chr\t$gene_start\t$gene_end\t$segment1_chr\t$segment1_start\t$segment1_end\t$segment2_chr\t$segment2_start\t$segment2_end";
        if(exists $hash1{$k}){
            my $right_overlap_bp =$hash1{$k};
            print $O1 "$k\t$left_overlap_bp\t$right_overlap_bp\n";
        }
    }
}

