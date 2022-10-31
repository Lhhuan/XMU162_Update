#"/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 注释ensg位置得"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot_0.05egene_pos.bed.gz",排序得"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos_sorted.bed.gz"

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my (%hash1,%hash2);

my $fo1 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos.bed.gz";
my $fo2 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos_sorted.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Hotspot_gene\t${type}\n";
my $f1= "/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2= "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
while(<$I1>)
{
    chomp;
    unless(/Gene/){
        my @f=split/\t/;
        my $ensg =$f[0];
        my $start =$f[1];
        my $end=$f[2];
        my $chr=$f[3];
        unless($chr=~/CHR|GL|KI|MT|X|Y/){
            my $v= "chr${chr}\t$start\t$end";
            $hash1{$ensg}=$v;
        }
    }
}

while(<$I2>)
{
    chomp;
    my @f=split/\t/;
    my $ensg=$f[3];
    if(exists $hash1{$ensg}){
        my $gene_pos = $hash1{$ensg};
        my $output ="$gene_pos\t$ensg";
        unless(exists $hash2{$output}){
            $hash2{$output}=1;
            print $O1 "$output\n";
        }

    }
    # else{
    #     print "$ensg\n";
    # }
}
close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > $fo2";