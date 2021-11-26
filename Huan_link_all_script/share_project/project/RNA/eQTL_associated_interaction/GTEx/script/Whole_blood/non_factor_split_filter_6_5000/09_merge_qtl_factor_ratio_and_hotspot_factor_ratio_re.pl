#---------------------------------------------------#对文件进行处理

# my @tissues = ("Whole_Blood","Lung");

my $tissue = "Whole_Blood";

$output_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/";
my $f1 = "$output_dir/03_${tissue}_count_eQTL_hit_factor.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/08_whole_Blood_prepare_number_ROC_factor_count_re.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my(%hash1,%hash2,%hash3,%hash4);


my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/09_Whole_Blood_merge_qtl_factor_ratio_and_hotspot_factor_ratio_re.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Factor/){
        my $Factor = $f[0];
        my $hit_ratio =$f[-1];
        $hash1{$Factor}=$hit_ratio;
    }
}


while(<$I2>){
    chomp;
    my @f =split/\t/;
    if(/^Factor/){
        print $O1 "$_\teQTL_hit_ratio\tsignificant\n";
    }
    else{
        my $Factor = $f[0];
        my $TPR =$f[-2];
        if (exists $hash1{$Factor}){
            my $eqtl_hit_ratio = $hash1{$Factor};
            if ($TPR >$eqtl_hit_ratio){
                my $sig = "T";
                print $O1 "$_\t$eqtl_hit_ratio\t$sig\n";
            }
            else{
                my $sig = "F";
                print $O1 "$_\t$eqtl_hit_ratio\t$sig\n";               
            }
        }
    }
}
