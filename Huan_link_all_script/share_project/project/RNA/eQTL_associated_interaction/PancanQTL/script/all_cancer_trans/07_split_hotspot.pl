#将hotspot makewindows 1MB  ../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}_makewin_${win}.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;

my $command2= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  261 | gzip > hg19.chrom1_22_sizes_sorted.261BP_windows.bed.gz";
system "$command2";
my @tissues =  ("BRCA","KIRC");
my $win = "261BP";
foreach my $tissue1(@tissues){
    my $tissue11 = $tissue1; 
    $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
    my $file1 = "../../output/${tissue1}/Trans_eQTL/hotspot_trans_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz"; 
    my $command2 = "bedtools intersect -a $file1 -b hg19.chrom1_22_sizes_sorted.${win}_windows.bed.gz -wo |cut -f4-7 |gzip > ../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}_makewin_${win}.bed.gz";
    system $command2;

}
