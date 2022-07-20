#产生10000个与../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/Tissue_merge_segment_${group}_cutoff_${cutoff}.bed.gz相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Tissue_merge";

# my $input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz"; #hotspot 
my $input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"; #hotspot 
my $input_file_base_name = basename($input_file);
# $input_file_base_name =~ s/Tissue_merge/Tissue_merge/g;
my $output_dir=  "../output";

my $whole_genome = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted.txt";

my $basic_simulate = "01_extend_5_fold_neighbour_sort_merge.bed.gz";

my $command10 = "bedtools complement -i ${output_dir}/$basic_simulate -g $whole_genome | sort -k1,1 -k2,2n |  gzip > ${output_dir}/complement_$basic_simulate"; #emp 0 的互补 
system $command10;

my $command11 = "zless ${output_dir}/complement_$basic_simulate > ${output_dir}/complement_final_5_fold_neighbour.bed ";

my $command12 = "zless $input_file >> ${output_dir}/complement_final_5_fold_neighbour.bed"; #Hotspot 进行了extend，所以exclude不仅仅是emp 0 的互补，还有extend的部分 既最终用不在hotspot中的emp0为背景

system $command11;
system $command12;
system "less ${output_dir}/complement_final_5_fold_neighbour.bed |sort -k1,1 -k2,2n |gzip >${output_dir}/complement_final_5_fold_neighbour_sorted.bed.gz ";
system "bedtools merge -i ${output_dir}/complement_final_5_fold_neighbour_sorted.bed.gz |gzip >${output_dir}/complement_final_5_fold_neighbour_sorted_merge.bed.gz"; # exclude  

my $genome="/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22.sizes";


# for (my $i=1;$i<10001;$i++){
    my $i=1;
    my $out_file = "$output_dir/${i}_resemble_${input_file_base_name}";
    #generate random file 
    my $command1 = "bedtools shuffle -i $input_file -g $genome -excl ${output_dir}/complement_final_5_fold_neighbour_sorted_merge.bed.gz -chrom | gzip >$out_file"; #即用$basic_simulate 除去hotspot extend 的部分 进行随机抽样
    # print "$command\n";
    system $command1;
#     print "$i\n";
# }



