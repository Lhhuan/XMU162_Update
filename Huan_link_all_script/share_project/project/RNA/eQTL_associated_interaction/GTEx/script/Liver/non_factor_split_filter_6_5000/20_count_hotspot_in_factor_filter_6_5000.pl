#根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/061_hotspot_in_different_cutoff_interval_18.txtt记录的hotspot数目，和annotation信息得每个cutoff下hit住factor的hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor.txt.gz"; 得每个cutoff下被每种factor hit住hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_per_factor_ratio.txt.gz";得每个cutoff下被每种factor hit住hotspot的文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor_pos.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
for (my $i=0.1;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}

my $tissue= "Lung";
my %hash1;
my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/061_hotspot_in_different_cutoff_interval_18.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
while(<$I1>)
{
    chomp;
    unless(/Cutoff/){
        my @f = split/\t/;
        my $Cutoff = $f[0];
        my $Number_of_hotspot = $f[1];
        $hash1{$Cutoff}=$Number_of_hotspot;
    }
}

my $type = "factor";
my $group = "hotspot";
my %hash9;


my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Cutoff\tnum_of_hotspot_hit_factor\tnum_all_hotspot\thit_ratio\n";
my $fo2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor_pos.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "chr\tstart\tend\tfactor\tcutoff\n";

my $fo3 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_per_factor_ratio.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;
print $O3 "cutoff\thot_hit_enhancer\thot_hit_enhancer_ratio\thot_hit_promoter\thot_hit_promoter_ratio\thot_hit_TFBS\thot_hit_TFBS_ratio\thot_hit_CHROMATIN_Accessibility\thot_hit_CHROMATIN_Accessibility_ratio\thot_hit_HISTONE_modification\thot_hit_HISTONE_modification_ratio\thot_hit_CTCF\thot_hit_CTCF_ratio\tall_unique_hotspot_hitted\tall_hotspot_num\tall_ratio\n";

foreach my $cutoff(@cutoffs){
    my $all_hotspot_num = $hash1{$cutoff};
    my $anno_dir = "../../../output/${tissue}/Cis_eQTL/annotation/interval_18_filter/6_5000/ALL/${type}/${group}/${cutoff}"; 
    my $base_file = "${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    my (%hash2,%hash3,%hash4,%hash5,%hash6,%hash7,%hash8);
    my $f2 = "${anno_dir}/enhancer_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash2{$pos}=1;
        push @{$hash8{$pos}},"enhancer";

    }

    my $f3 = "${anno_dir}/promoter_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
    while(<$I3>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash3{$pos}=1;
        push @{$hash8{$pos}},"promoter";
    }

    my $f4 = "${anno_dir}/TFBS_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
    while(<$I4>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash4{$pos}=1;
        push @{$hash8{$pos}},"TFBS";
    }

    my $f5 = "${anno_dir}/CHROMATIN_Accessibility_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件
    while(<$I5>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash5{$pos}=1;
        push @{$hash8{$pos}},"CHROMATIN_Accessibility";
    }

    my $f6 = "${anno_dir}/HISTONE_modification_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I6 ,"gzip -dc $f6|") or die ("can not open input file '$f6' \n"); #读压缩文件
    while(<$I6>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash6{$pos}=1;
        push @{$hash8{$pos}},"HISTONE_modification";
    }

    my $f7 = "${anno_dir}/CTCF_${base_file}";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I7 ,"gzip -dc $f7|") or die ("can not open input file '$f7' \n"); #读压缩文件
    while(<$I7>)
    {
        chomp;
        my @f = split/\t/;
        my $pos = join("\t",@f[0..2]);
        $hash7{$pos}=1;
        push @{$hash8{$pos}},"CTCF";
    }

    my $num_factor_hotspot = keys %hash8;
    my $ratio = $num_factor_hotspot/$all_hotspot_num;
    my $output = "$cutoff\t$num_factor_hotspot\t$all_hotspot_num\t$ratio";
    # print $O1 "$cutoff\t$num_factor_hotspot\t$all_hotspot_num\t$ratio\n";
    unless(exists $hash9{$output}){
        $hash9{$output} =1;
        print $O1 "$output\n";
    }

    foreach my $pos (sort keys %hash8){
        my @factors = @{$hash8{$pos}};
        my %hash;
        @factors = grep { ++$hash{$_} < 2 } @factors;
        my $factor_hotspot = join(";",@factors);

        print $O2 "$pos\t$factor_hotspot\t$cutoff\n";
    }
    my $enhancer_h = keys %hash2; #获得哈希长度
    my $promoter_h = keys %hash3;
    my $TFBS_h = keys %hash4;
    my $CHROMATIN_Accessibility_h = keys %hash5;
    my $HISTONE_modification_h = keys %hash6;
    my $CTCF_h = keys %hash7;
    my $enhancer_h_ratio = $enhancer_h/$all_hotspot_num;
    my $promoter_h_ratio = $promoter_h/$all_hotspot_num;
    my $TFBS_h_ratio = $TFBS_h/$all_hotspot_num;
    my $CHROMATIN_Accessibility_h_ratio = $CHROMATIN_Accessibility_h/$all_hotspot_num;
    my $HISTONE_modification_h_ratio = $HISTONE_modification_h/$all_hotspot_num;
    my $CTCF_h_ratio = $CTCF_h/$all_hotspot_num;
    print $O3 "$cutoff\t$enhancer_h\t$enhancer_h_ratio\t$promoter_h\t$promoter_h_ratio\t$TFBS_h\t$TFBS_h_ratio\t$CHROMATIN_Accessibility_h\t$CHROMATIN_Accessibility_h_ratio\t$HISTONE_modification_h\t$HISTONE_modification_h_ratio\t$CTCF_h\t$CTCF_h_ratio\t$num_factor_hotspot\t$all_hotspot_num\t\t$ratio\n";
}

