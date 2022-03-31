#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;




my $genome="/share/Projects/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
my $chr1_hotspot = "/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/script/chr1/all/output/07_all_predict_need.bed.gz";
my $exclude_file = "/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz";
my $out_file= "../output/07_random_chr1.bed.gz";
# my $command1 = "bedtools shuffle -i $chr1_hotspot -g $genome -excl $chr1_hotspot -chrom | gzip >$out_file"; #即用$basic_simulate 除去hotspot extend 的部分 进行随机抽样
my $command1 = "bedtools shuffle -i $chr1_hotspot -g $genome -excl $exclude_file -chrom | gzip >$out_file"; #即用$basic_simulate 除去hotspot extend 的部分 进行随机抽样
system $command1;