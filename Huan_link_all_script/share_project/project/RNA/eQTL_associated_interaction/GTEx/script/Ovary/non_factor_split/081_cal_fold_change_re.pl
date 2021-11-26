#---------------------------------------------------#对文件进行处理

# my @tissues = ("Ovary","Ovary");

# my $tissue = "Ovary";

# $output_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/";
my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Ovary/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_Ovary_prepare_number_ROC_factor_count_re.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

# my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Ovary/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_Ovary_prepare_number_ROC_factor_count.txt";
# # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
# my(%hash1,%hash2,%hash3,%hash4);


my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Ovary/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_Ovary_cal_fold_change_re.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Factor\tCutoff\tNumber_of_factor_in_hotspot_TP\tNumber_of_factor_in_non_hotspot_FN\tFold_change_value\n";

while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Factor/){
        my $Factor = $f[0];
        my $Cutoff =$f[1];
        my $TP = $f[2];
        my $FN = $f[3];
        my $fold_change = $TP/$FN;
        print $O1 "$Factor\t$Cutoff\t$TP\t$FN\t$fold_change\n";
    }
}
