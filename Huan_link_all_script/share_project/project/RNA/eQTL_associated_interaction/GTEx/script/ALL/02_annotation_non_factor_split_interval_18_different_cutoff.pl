#对"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/significant_eQTL.txt.gz"进行anno得../../output/${tissue}/Cis_eQTL/eQTL_annotation/
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my @tissues = ("Whole_Blood","Lung");
my @tissues = ("Whole_Blood","Lung","Ovary","Liver","Pancreas");

# my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
my $type = "factor";
foreach my $tissue(@tissues){
    # my $pid = $pm->start and next; #开始多线程
    my $input_file = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/significant_eQTL.txt.gz";
    # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
    my $input_file_base_name = basename($input_file);
    # my $dir = dirname($script);
    # print "$input_file_base_name\n";
    my $output_dir = "../../output/${tissue}/Cis_eQTL/eQTL_annotation/"; 
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
    $ENV{'tissue'} = $tissue ;
    my $command = "bash annotation_${type}_bedtools_intersect_interval18.sh";
    system $command;
    print "$tissue\n";
    # $pm->finish;  #多线程结束
}

