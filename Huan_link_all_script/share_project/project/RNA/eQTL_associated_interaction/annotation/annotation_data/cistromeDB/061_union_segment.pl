#将normal cell ${dir}/merge_pos_info_narrow_peak_sort.bed.gz bedtools merge -i 为${dir}/merge_pos_info_narrow_peak_sort_union.bed.gz，得互补文件${dir}/merge_pos_info_narrow_peak_union_complement.bed.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;



# my @files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt");
my @dirs = ("./normal_cell/Human_FACTOR","./normal_cell/HISTONE_MARK_AND_VARIANT","./normal_cell/Human_CHROMATIN_Accessibility");
# my $pm = Parallel::ForkManager->new(4); ## 设置最大的线程数目

foreach my $dir(@dirs){
    my $input = "${dir}/merge_pos_info_narrow_peak_sort.bed.gz";
    my $output= "${dir}/merge_pos_info_narrow_peak_sort_union.bed.gz";
    my $output2= "${dir}/merge_pos_info_narrow_peak_union_complement.bed.gz";
    my $command = "bedtools merge -i $input |gzip > $output";
    my $command2 = "bedtools  complement -i $output -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt | gzip> $output2";
    # print "$command\n";
    system $command;
    system $command2;
    print "$dir\n";
}