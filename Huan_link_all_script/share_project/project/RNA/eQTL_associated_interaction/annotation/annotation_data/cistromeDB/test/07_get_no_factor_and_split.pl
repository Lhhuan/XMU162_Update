#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;




my @dirs = ("./normal_cell/Human_FACTOR","./normal_cell/HISTONE_MARK_AND_VARIANT","./normal_cell/Human_CHROMATIN_Accessibility");
foreach my $dir(@dirs){
    my $input = "$dir/merge_pos_info_narrow_peak.bed.gz";
    my $input_sort = "$dir/merge_pos_info_narrow_peak_sort.bed.gz";
    my $command1 = "zless $input |sort -k1,1 -k2,2n |gzip >$input_sort";
    system $command1;
    my $output = "$dir/merge_pos_info_narrow_peak_complement.bed.gz";
    my $command2 = "bedtools  complement -i $input_sort -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt | gzip> $output";
    # print "$command\n";
    
    system $command2;
}