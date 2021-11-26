#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/11_F-score.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Factor\tCutoff\tNumber_of_factor_in_hotspot_TP\tNumber_of_factor_in_non_hotspot_FN\tNumber_of_non_factor_in_hotspot_FP\tNumber_of_non_factor_in_non_hotspot_TN\tTPR\tFPR\tPrecision\trecall\tF_score\n";



while(<$I1>)
{
    chomp;
    unless(/Factor/){
        my @f = split/\t/;
        my $Factor = $f[0];
        my $Cutoff = $f[1];
        if($Cutoff <=0.65){
            my $TP= $f[2];
            my $FN = $f[3];
            my $FP = $f[4];
            my $TN = $f[5];
            # my $v= "$pos\t$emplambda";
            my $precision = $TP/($TP+$FP);
            my $recall = $TP/($TP+$FN);
            my $F_score = 2*($precision*$recall)/($precision + $recall);
            $F_score = sprintf "%.5f",$F_score;
            print $O1 "$_\t$precision\t$recall\t$F_score\n";
        }
    }
}