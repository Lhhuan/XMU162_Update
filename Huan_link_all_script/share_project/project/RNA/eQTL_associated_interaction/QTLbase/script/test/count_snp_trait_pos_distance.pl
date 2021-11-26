# 计算 ../output/01_all_kinds_QTL.txt中snp 与trait的位置差，得../output/count_snp_trait_pos_distance.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../output/01_all_kinds_QTL.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../output/count_snp_trait_pos_distance.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/count_snp_trait_pos_in_diff_chr.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $header = "SNP_chr\tSNP_pos\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tPvalue\tSourceid\tQTL_type\tSNP-Trait_distance\tabs_SNP-Trait_distance";
print  $O1 "$header\n";
print  $O2 "SNP_chr\tSNP_pos\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tPvalue\tSourceid\tQTL_type\n";
#-------------------------------------------------------------#获取组织名称


while(<$I1>)
{
    chomp;
    unless(/SNP_chr/){
        my @f =split/\t/;
        my $SNP_chr =$f[0];
        my $SNP_pos =$f[1];
        my $Trait_chr =$f[3];
        my $Trait_start =$f[4];
        my $Trait_end =$f[5];
        if ($SNP_chr == $Trait_chr){
            # print "$SNP_chr\t$Trait_chr\n";
            my $Trait_median = ($Trait_start + $Trait_end)/2;
            $Trait_median =int($Trait_median);
            my $distance = $Trait_median - $SNP_pos;
            my $abs_distance = abs($distance);
            print $O1 "$_\t$distance\t$abs_distance\n";
        }
        else{ #不在同一染色体上
            print $O2 "$_\tDiff_chr\n";
        }
    }
}
