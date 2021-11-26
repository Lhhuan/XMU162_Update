#根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/annotation/interval_18/ALL/${number}/factor/${group}/${cutoff}/${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"，"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/annotation/interval_18/ALL/${number}/non_factor_split/${group}/${cutoff}/non_${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz" 以factor为基数，准备计算ROC的四格表得""/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/${tissue}/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
# my @types = ("factor","non_factor");
my @factors = ("promoter","enhancer","TFBS","CHROMATIN_Accessibility","HISTONE_modification","CTCF");
my @tissues = ("Whole_Blood","Lung","Ovary","Liver","Pancreas");


foreach my $tissue(@tissues){
    my $output_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }

    my $fo1 = "$output_dir/03_${tissue}_count_eQTL_hit_factor.txt";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

    print $O1 "Factor\teQTL_hit_factor_count\tfactor_count\thit_ratio\n";
    # TP\t$FN\t$FP\t$TN

    foreach my $factor(@factors){
        my $command_factor = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/${factor}_union.bed.gz | wc -l" ;
        my $factor_line_count = wc($command_factor);
        my $hit_factor_file = "../../output/${tissue}/Cis_eQTL/eQTL_annotation/${factor}_significant_eQTL.txt.gz";
        #--------------------------------------get size of file 
        my @arg1s = stat ($hit_factor_file);
        my $factor_size = $arg1s[7];
        #-----------------------------
        if($factor_size >20){
             my $command_factor_hit = "zless $hit_factor_file |cut -f4,5,6| sort -u | wc -l" ;
             my $factor_hit_count = wc($command_factor_hit);
             my $hit_ratio = $factor_hit_count/$factor_line_count;
             print $O1 "$factor\t$factor_hit_count\t$factor_line_count\t$hit_ratio\n";
        }
    }
}

sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
