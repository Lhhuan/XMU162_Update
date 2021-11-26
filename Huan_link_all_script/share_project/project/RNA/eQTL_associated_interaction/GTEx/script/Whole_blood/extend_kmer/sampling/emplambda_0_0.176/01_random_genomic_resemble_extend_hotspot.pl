#产生1000个与"../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6/extend_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"相同的resemble hotspot,"../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random_extend/6/background_original_random/${cutoff}/$cutoff2"
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
my $tissue = "Whole_Blood";

my $input_file = "../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6/extend_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
my $input_file_base_name = basename($input_file);
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
my $cutoff2 = "0_0.176";
my $output_dir=  "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/$cutoff2/bed";

if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}

# my $cutoff2 = "0_0.176";
my $anno_dir = "../../../../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/sampling";
my $whole_genome = "/home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt";

my $basic_simulate = "whole_blood_segment_hotspot_cutoff_0.bed.gz";

my $command10 = "bedtools complement -i $anno_dir/$basic_simulate -g $whole_genome | sort -k1,1 -k2,2n |  gzip > $anno_dir/complement_$basic_simulate"; #exclude
system $command10;


my $genome="/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";


for (my $i=1;$i<1001;$i++){
    # my $i=1;
    my $out_file = "$output_dir/${i}_resemble_${input_file_base_name}";
    #generate random file 
    my $command1 = "bedtools shuffle -i $input_file -g $genome -excl $anno_dir/complement_$basic_simulate -chrom | gzip >$out_file"; #即用$basic_simulate 进行随机抽样
    # print "$command\n";
    system $command1;
    my $fa_out_dir = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/$cutoff2/fa";
    my $out_file_base_name = "${i}_resemble_${input_file_base_name}";
    $out_file_base_name =~ s/bed\.gz/fa/;
    my $command2 = " bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed  $out_file -fo $fa_out_dir/$out_file_base_name";
    system $command2;
    print "$i\n";
}



