#@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时,，对interval为15_6_7_8_9_12_18_时，用041_bedtools_intersect_hotspot_overlap.sh进行overlap进行overlap,
#得../../../output/ALL_eQTL/cis_trans/$group/segment_overlap/15_6_7_8_9_12_18_cutoff_7.3_${type}_eQTL_segment_$group.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $QTL = "eQTL";
my @types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my @groups = ("hotspot","non_hotspot");
my $pm = Parallel::ForkManager->new(2); ## 设置最大的线程数目
foreach my $type(@types){
    my $pid = $pm->start and next; #开始多线程
    foreach my $group(@groups){
        print "$group\t$type\n";

        $ENV{'type'}  = $type; #设置环境变量
        $ENV{'group'} = $group ;
        my $command = "bash 041_bedtools_intersect_hotspot_overlap.sh";
        # print "$command\n";
        system $command;
    } 
    $pm->finish;  #多线程结束   
}


