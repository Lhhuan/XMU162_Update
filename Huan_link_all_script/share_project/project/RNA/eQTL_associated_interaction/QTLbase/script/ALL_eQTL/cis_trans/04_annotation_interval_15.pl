#@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时,用annotation_bedtools_intersect.sh进行annotation,得$output_dir/RBP_$input_file_base_name
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $QTL = "eQTL";
my $interval = 15;
my @types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my @groups = ("hotspot","non_hotspot");
my $pm = Parallel::ForkManager->new(2); ## 设置最大的线程数目
foreach my $type(@types){
    my $pid = $pm->start and next; #开始多线程
    foreach my $group(@groups){
        my $input_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
        # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
        my $input_file_base_name = basename($input_file);
        # my $dir = dirname($script);
        # print "$input_file_base_name\n";
        my $output_dir = "../../../output/ALL_eQTL/cis_trans/interval_15/annotation/${group}/${type}/";

        $ENV{'input_file'}  = $input_file; #设置环境变量
        $ENV{'input_file_base_name'} = $input_file_base_name ;
        $ENV{'output_dir'} = $output_dir ;
        my $command = "bash annotation_bedtools_intersect.sh";
        # print "$command\n";
        system $command;
    } 
    $pm->finish;  #多线程结束   
}


