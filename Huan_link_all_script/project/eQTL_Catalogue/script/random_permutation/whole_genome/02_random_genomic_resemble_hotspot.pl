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

my $input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz"; #hotspot 
my $input_file_base_name = basename($input_file);
# $input_file_base_name =~ s/Tissue_merge/Tissue_merge/g;
my $output_dir=  "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend/whole_genome";

if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}

my $genome="/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22.sizes";


for (my $i=1;$i<10001;$i++){
    # my $i=1;
    my $out_file = "$output_dir/${i}_resemble_${input_file_base_name}";
    #generate random file 
    my $command1 = "bedtools shuffle -i $input_file -g $genome -excl $input_file -chrom | gzip >$out_file"; #即用$basic_simulate 除去hotspot extend 的部分 进行随机抽样
    # print "$command\n";
    system $command1;
    print "$i\n";
}


