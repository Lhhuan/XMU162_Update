#将#/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor_summary.txt和
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt 
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt merge 在一起
#得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/12_merge_random_and_original_fisher_exact_test.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/12_merge_random_and_original_fisher_exact_test.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

print $O1 "Tissue\tFactor\tCutoff\t_ran_sum_Pvalue\tran_sum_OR\tran_sum_Significant\ttrue_Pvalue\ttrue_OR\ttrue_Significant\n";


my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

my $f3 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor_summary.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 

my(%hash1,%hash2,%hash3,%hash4);


#----------------------
while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Pvalue/){
        my $Pvalue = $f[0];
        my $OR =$f[1];
        my $Significant =$f[2];
        my $Factor =$f[3];
        my $Cutoff =$f[4];
        my $k = "Lung\t$Factor\t$Cutoff";
        my $v = join ("\t",@f[0..2]);
        $hash1{$k}=$v;
    }
}


#-------------------

while(<$I2>){
    chomp;
    my @f = split/\t/;
    unless(/^Pvalue/){
        my $Pvalue = $f[0];
        my $OR =$f[1];
        my $Significant =$f[2];
        my $Factor =$f[3];
        my $Cutoff =$f[4];
        my $k = "Whole_Blood\t$Factor\t$Cutoff";
        my $v = join ("\t",@f[0..2]);
        $hash1{$k}=$v;
    }
}



while(<$I3>){
    chomp;
    my @f = split/\t/;
    unless(/^Tissue/){
        my $Tissue = $f[0];
        my $Random_number = $f[1];
        my $Pvalue = $f[2];
        my $OR =$f[3];
        my $Significant =$f[4];
        my $Factor =$f[5];
        my $Cutoff =$f[6];
        my $k = "$Tissue\t$Factor\t$Cutoff";
        if (exists $hash1{$k}){
            my $v= $hash1{$k};
            print $O1 "$k\t$Pvalue\t$OR\t$Significant\t$v\n";
        }
    }
}
