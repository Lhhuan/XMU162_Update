#将/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt中random的数据合在一起，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

print $O1 "Tissue\tRandom_number\tFactor\tCutoff\tNumber_of_factor_in_hotspot_TP\tNumber_of_factor_in_non_hotspot_FN\tNumber_of_non_factor_in_hotspot_FP\tNumber_of_non_factor_in_non_hotspot_TN\tTPR\tFPR\n";


my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my(%hash1,%hash2,%hash3,%hash4);
while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Tissue/){
        my $Tissue = $f[0];
        my $Random_number = $f[1];
        my $Factor = $f[2];
        my $Cutoff =$f[3];
        my $Number_of_factor_in_hotspot_TP =$f[4];
        my $Number_of_factor_in_non_hotspot_FN =$f[5];
        my $Number_of_non_factor_in_hotspot_FP =$f[6];
        my $Number_of_non_factor_in_non_hotspot_TN =$f[7];
        my $k = "$Tissue\t$Factor\t$Cutoff";
        push @{$hash1{$k}},$Number_of_factor_in_hotspot_TP;
        push @{$hash2{$k}},$Number_of_factor_in_non_hotspot_FN;
        push @{$hash3{$k}},$Number_of_non_factor_in_hotspot_FP;
        push @{$hash4{$k}},$Number_of_non_factor_in_non_hotspot_TN;
    }
}


foreach my $k(sort keys %hash1){
    my @TPs = @{$hash1{$k}};
    my @FNs = @{$hash2{$k}};
    my @FPs = @{$hash3{$k}};
    my @TNs = @{$hash4{$k}};
    my $TP = sum @TPs;
    my $FN = sum @FNs;
    my $FP = sum @FPs;
    my $TN = sum @TNs;
    my $TPR =$TP/($TP+$FN);
    my $FPR = $FP/($FP +$TN);
    my @t =split/\t/,$k;
    my $Tissue =$t[0];
    my $Factor =$t[1];
    my $Cutoff =$t[2];
    my $number = "all";
    print $O1 "$Tissue\t$number\t$Factor\t$Cutoff\t$TP\t$FN\t$FP\t$TN\t$TPR\t$FPR\n";
    
}