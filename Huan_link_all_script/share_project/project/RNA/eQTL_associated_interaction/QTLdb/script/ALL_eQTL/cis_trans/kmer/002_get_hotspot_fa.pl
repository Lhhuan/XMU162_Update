#用 ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa 为"../../../../output/ALL_eQTL/cis_trans/hotspot/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.bed.gz" 提取序列得 ,用a"${output_dir}/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot.fa"
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
my @intervals = (12,15,18);
my @types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my $group = "hotspot";
my $pm = Parallel::ForkManager->new(6); ## 设置最大的线程数目
foreach my $interval(@intervals){
    my $pid = $pm->start and next; #开始多线程
    foreach my $type(@types){
        my $output_dir = "../../../../output/ALL_eQTL/cis_trans/kmer/interval_${interval}/${type}";
        if(-e $output_dir){
            print "${output_dir}\texist\n";
        }
        else{
            system "mkdir -p $output_dir";
        }
        my $input_file = "../../../../output/ALL_eQTL/cis_trans/hotspot/random_select/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.bed.gz";
        my $output_file = "${output_dir}/interval_${interval}_cutoff_7.3_${type}_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.fa";
        my $command = "bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed ${input_file} -fo $output_file";
        # print "$command\n";
        system $command;
    }
    $pm->finish;  #多线程结束
}

