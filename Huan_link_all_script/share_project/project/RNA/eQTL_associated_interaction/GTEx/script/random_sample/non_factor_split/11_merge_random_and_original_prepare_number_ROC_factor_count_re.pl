#/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result.txt和
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt 
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt merge 在一起
#得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/11_merge_random_and_original_prepare_number_ROC_factor_count.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/11_merge_random_and_original_prepare_number_ROC_factor_count_re.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

print $O1 "Tissue\tRandom_number\tFactor\tCutoff\trandom_sum_Number_of_factor_in_hotspot_TP\trandom_sum_Number_of_factor_in_non_hotspot_FN\trandom_sum_Number_of_non_factor_in_hotspot_FP\trandom_sum_Number_of_non_factor_in_non_hotspot_TN\trandom_sum_TPR\trandom_sum_FPR\ttrue_Number_of_factor_in_hotspot_TP\ttrue_Number_of_factor_in_non_hotspot_FN\ttrue_Number_of_non_factor_in_hotspot_FP\ttrue_Number_of_non_factor_in_non_hotspot_TN\ttrue_TPR\ttrue_FPR\tTPR_percentage_true_greater_than_random\n";


my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_Lung_prepare_number_ROC_factor_count_re.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_whole_Blood_prepare_number_ROC_factor_count_re.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

my $f3 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result_re.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 

my(%hash1,%hash2,%hash3,%hash4);


#----------------------
while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Factor/){
        my $Factor = $f[0];
        my $Cutoff =$f[1];
        my $Number_of_factor_in_hotspot_TP =$f[2];
        my $Number_of_factor_in_non_hotspot_FN =$f[3];
        my $Number_of_non_factor_in_hotspot_FP =$f[4];
        my $Number_of_non_factor_in_non_hotspot_TN =$f[5];
        my $TPR = $f[-2];
        my $FPR = $f[-1];
        my $k = "Lung\t$Factor\t$Cutoff";
        my $v = join ("\t",@f[2..$#f]);
        $hash1{$k}=$v;
    }
}


#-------------------
while(<$I2>){
    chomp;
    my @f = split/\t/;
    unless(/^Factor/){
        my $Factor = $f[0];
        my $Cutoff =$f[1];
        my $Number_of_factor_in_hotspot_TP =$f[2];
        my $Number_of_factor_in_non_hotspot_FN =$f[3];
        my $Number_of_non_factor_in_hotspot_FP =$f[4];
        my $Number_of_non_factor_in_non_hotspot_TN =$f[5];
        my $TPR = $f[-2];
        my $FPR = $f[-1];
        my $k = "Whole_Blood\t$Factor\t$Cutoff";
        my $v = join ("\t",@f[2..$#f]);
        $hash1{$k}=$v;
    }
}
#--------------



while(<$I3>){
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
        my $TPR = $f[-2];
        my $k = "$Tissue\t$Factor\t$Cutoff";
        if (exists $hash1{$k}){
            my $v= $hash1{$k};
            my @t = split/\t/,$v;
            my $true_TPR = $t[-2];
            my $add_ratio = ($true_TPR-$TPR)/$true_TPR *100;
            my $var = sprintf "%.2f",$add_ratio;
            print $O1 "$_\t$v\t${var}%\n";
        }
    }
}
