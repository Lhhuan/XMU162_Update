#interval_18 时，"../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz"用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
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
my $type = "factor";
my $tissue ="Whole_Blood";

my $input_file = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz";

my $input_file_base_name = basename($input_file);
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
# my $dir = dirname($script);
# print "$input_file_base_name\n";
my $output_dir = "../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${type}/${group}/${cutoff}"; 
# mkdir $PMID;
# #------------
if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}
#------------
$ENV{'input_file'}  = $input_file; #设置环境变量
$ENV{'input_file_base_name'} = $input_file_base_name ;
$ENV{'output_dir'} = $output_dir ;
# $ENV{'fraction'} = $fraction ;
my $command = "bash annotation_${type}_bedtools_intersect_interval18.sh";
# print "$command\n";
system $command;
print "$type\t$group\t$cutoff\n";
          
