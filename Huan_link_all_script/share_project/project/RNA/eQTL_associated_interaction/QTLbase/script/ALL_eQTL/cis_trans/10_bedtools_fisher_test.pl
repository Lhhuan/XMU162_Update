#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
# use Parallel::ForkManager;

# my $fraction = 0.5;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/fisher_exact_test/fisher_exact_test.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "In_a_in_b\tIn_a_not_in_b\tNot_in_a_in_b\tNot_in_a_not_in_b\tP_left\tP_right\tP_two_tail\tOdd_ratio\tInterval\tOverlap_fraction\tType\tFactor\n";
my @fractions = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1);
my @factors = ("circRNA/03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz","cistromeDB/HISTONE_MARK_AND_VARIANT/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz","cistromeDB/Human_CHROMATIN_Accessibility/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz",
    "cistromeDB/Human_FACTOR/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz", "CTCF/03_all_unique_merge_ctcf_binding_site_narrow_peak_sorted_chr1_22.bed.gz","enhancer/03_fantom5_enhancers_phase1_phase2_sorted_chr1_22.bed.gz",
    "hic/03_all_unique_hic_loops_tran_cis_1_2_sorted_chr1_22.bed.gz","lncRNA/03_NONCODEv5_hg19.lncAndGene_sorted_chr1_22.bed.gz","promoter/03_fantom5_promoter_phase1_phase2_sorted_chr1_22.bed.gz","RBP/03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz");
my @types = ("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my @interval = (6,7,8,9,12,15,18);
foreach my $j (@interval){
    foreach my $fraction(@fractions){
        foreach my $type(@types){
            foreach my $factor(@factors){
                print "$fraction\t$type\t$j\t$factor\n";
                # print "$fraction\n";
                my $command = "bedtools fisher -f $fraction  -a ../../../annotation_data/$factor -b ../../../output/ALL_eQTL/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_eQTL_segment_hotspot_length_of_segment.bed.gz   -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt";
                # my $command2 = "bedtools fisher -a ../../../annotation_data/circRNA/03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz -b ../../../output/ALL_eQTL/cis_trans/hotspot/interval_12_cutoff_7.3_cis_10MB_eQTL_segment_hotspot_length_of_segment.bed.gz   -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt";
                # print "$command\n";
                # system "$command";
                # system "$command2";
                my $results = readpipe($command);
                my @f = split/\n/,$results;
                my $aa = $f[8];
                $aa =~ s/\s+//g;
                $aa =~ s/#in-a|//g; 
                my @aass = split/\|/,$aa;
                my $in_a_in_b = $aass[1];
                my $in_a_not_in_b =$aass[2];
                my $bb = $f[9];
                $bb =~ s/\s+//g;
                $bb =~ s/#in-a|//g;
                my @bbss = split/\|/,$bb;
                my $not_in_a_in_b = $bbss[1];
                my $not_in_a_not_in_b = $bbss[2];
                my $P = $f[13];
                my @Ps = split/\s+/,$P;
                my $P_left = $Ps[0];
                my $P_right = $Ps[1];
                my $P_two_tail = $Ps[2];
                my $ratio = $Ps[3];
                my $output1 =  "$in_a_in_b\t$in_a_not_in_b\t$not_in_a_in_b\t$not_in_a_not_in_b";
                my $output2 = "$P_left\t$P_right\t$P_two_tail\t$ratio";
                # my @fff = split/\//,$factor;
                # print "@fff\n";
                $factor =~ s/\/03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz|\/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz|\/03_all_unique_merge_ctcf_binding_site_narrow_peak_sorted_chr1_22.bed.gz|\/03_fantom5_enhancers_phase1_phase2_sorted_chr1_22.bed.gz//g;
                $factor =~ s/\/03_all_unique_hic_loops_tran_cis_1_2_sorted_chr1_22.bed.gz|\/03_NONCODEv5_hg19.lncAndGene_sorted_chr1_22.bed.gz|\/03_fantom5_promoter_phase1_phase2_sorted_chr1_22.bed.gz|\/03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22.bed.gz//g;
                $factor =~ s/cistromeDB\///g;
                print "$factor\n";

                print $O1 "$output1\t$output2\t$j\t$fraction\t$type\t$factor\n";
            }
        }
    }
}