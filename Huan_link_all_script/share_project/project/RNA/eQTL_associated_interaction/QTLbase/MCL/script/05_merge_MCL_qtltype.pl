#将../output/04_part_mcl_result.txt.gz 和../../../output/01_all_kinds_QTL.txt.gz中特定QTL取overlap得 ../output/05_merge_MCL_qtltype_result.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
use File::Basename;

my $f1 = "../output/04_part_mcl_result.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../../output/01_all_kinds_QTL.txt.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/05_merge_MCL_qtltype_result.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "SNP\tBand\tCluster_number\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tQTL_type\n";
my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/^Cluster_number/){
        my @f = split/\t/;
        my $Cluster_number =$f[0];
        my $snp = $f[1];
        my $band =$f[2];
        my $v= "$band\t$Cluster_number";
        $hash1{$snp}=$v;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^SNP/){
        my @f = split/\t/;
        my $SNP_chr = $f[0];
        my $SNP_pos=$f[1];
        my $Mapped_gene =$f[2];
        my $Trait_chr =$f[3];
        my $Trait_start =$f[4];
        my $Trait_end =$f[5];
        my $Pvalue =$f[-3];
        my $QTL_type =$f[-1];
        if($QTL_type =~/\beQTL|mQTL|miQTL\b/){
            if ($Pvalue < 5e-8){
                my $k = "${SNP_chr}_${SNP_pos}";
                if (exists $hash1{$k}){
                    my $v = $hash1{$k};
                    print $O1 "$k\t$v\t$Mapped_gene\t$Trait_chr\t$Trait_start\t$Trait_end\t$QTL_type\n";
                }
            }
        }
    }
}

