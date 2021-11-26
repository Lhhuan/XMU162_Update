#@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot"),@intervals = (6,7,8,9,12,15,18),@fractions = (0.5,0.6,0.7,0.8,0.9,1)时,用annotation_bedtools_intersect_v2.sh进行annotation,得$output_dir/RBP_$input_file_base_name
#extend 100bp
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $QTL = "eQTL";
# my $interval = 18;
# my @intervals = (6,7,8,9,12,15,18);
my @intervals = (18,15);
# my @intervals = (18,15,12);
my @types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my @groups = ("hotspot","non_hotspot");
my @fractions = (0.1,1);
# my @fractions = (0.1);
my $pm = Parallel::ForkManager->new(2); ## 设置最大的线程数目
foreach my $interval(@intervals){
    foreach my $type(@types){
        foreach my $group(@groups){
            foreach my $fraction(@fractions){
                my $pid = $pm->start and next; #开始多线程
                my $input_file = "../../../output/ALL_${QTL}/cis_trans/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
                # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
                my $input_file_base_name = basename($input_file);
                # my $dir = dirname($script);
                # print "$input_file_base_name\n";
                my $output_dir = "/state/partition1/huan/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test_extend_500bp/annotation/interval_${interval}/${group}/${type}/${fraction}"; #compute-0-1
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
                $ENV{'fraction'} = $fraction ;
                my $command = "bash annotation_bedtools_intersect_v3.sh";
                # print "$command\n";
                system $command;
                print "$interval\t$type\t$group\t$fraction\n";
                $pm->finish;  #多线程结束
            }
        } 
    }
}

