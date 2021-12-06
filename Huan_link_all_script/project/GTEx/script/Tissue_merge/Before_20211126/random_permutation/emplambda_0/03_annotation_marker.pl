#对10000个../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2/*进行histone marker的注释，得/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2/*
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
# foreach my $type(@types){
my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Tissue_merge";
my $cutoff2 = "0";
my $output_dir =  "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/interval_18_random/original_random/emp${cutoff}/${cutoff2}";

my $anno_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/";
my $input_dir = $output_dir;
for (my $i=1;$i<10001;$i++){
# for (my $i=500;$i<1001;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $input_file_base_name = "${i}_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_sorted.bed.gz"; 
    my $input_file = "$input_dir/$input_file_base_name";
    my $sorted_input_file = "$input_dir/sorted_$input_file_base_name";
    my $command1 = "zless $input_file | sort -k1,1 -k2,2n | gzip > $sorted_input_file";
    $ENV{'sorted_input_file'}  = $sorted_input_file; #设置环境变量
    $ENV{'input_file_base_name'} = $input_file_base_name ;
    $ENV{'output_dir'} = $output_dir ;
    $ENV{'anno_dir'} = $anno_dir ;
    # my $out_basename = 
    # print "$i\n";
    # #-----------annotation
    my $command2 = "bash annotation_marker_interval18.sh";
    # my $command4 = "rm $input_file ";
    system $command1;
    system $command2;
    # system $command3;
    # system $command4;

    print "$i\n";
    $pm->finish;  #多线程结束
}


