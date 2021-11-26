#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
# use Parallel::ForkManager;

my $fraction = 0.9;
print "$fraction\n";
# my $command = "bedtools fisher -f $fraction  -a ../../../annotation_data/circRNA/03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz -b ../../../output/ALL_eQTL/cis_trans/hotspot/interval_12_cutoff_7.3_cis_10MB_eQTL_segment_hotspot_length_of_segment.bed.gz   -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt";
my $command2 = "bedtools fisher -f $fraction -a ../../../annotation_data/RBP/03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz -b ../../../output/ALL_eQTL/cis_trans/hotspot/interval_18_cutoff_7.3_trans_10MB_eQTL_segment_hotspot_length_of_segment.bed.gz   -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt";
# print "$command\n";
# system "$command";
system "$command2";
# my $results = readpipe($command);
# my @f = split/\n/,$results;
# print "name\n";
# print "$f[9]\n";
# print "$results\n";


# ../../../annotation_data/CTCF/03_all_unique_merge_ctcf_binding_site_narrow_peak_sorted_chr1_22.bed.gz

