#将用 ../data/cytoBand.txt.gz对compute-0-6下的/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz 
#进行annotation,得../output/03_annotation_ld_band.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my $f1 = "../data/cytoBand.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/03_annotation_ld_band.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/03_annotation_ld_not_band_info.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
my $header ="SNP_A\tSNP_B\tR2";
print $O1 "$header\tbind\n";
print $O2 "$header\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];
    $chr =~ s/chr//g;
    my $start =$f[1];
    my $end=$f[2];
    my $band =$f[3];
    my $k = "$chr\t$start\t$end";
    unless($chr=~/X|Y/){
        $hash1{$k}=$band;
        # print "$chr\n";
    }
    
}



while(<$I2>)
{
    chomp;
    unless (/CHR_A/){
        my $info =$_;
        $info =~s/^\s+//g;
        # print "$info\n";
        my @f = split/\s+/,$info;
        my $chr1 = $f[0];
        my $snp1 = $f[1];
        my $chr2 = $f[3]; #已确定$chr1 = $chr2
        my $snp2 = $f[4];
        my $snp_1 =$f[2];
        my $snp_2 = $f[-2];
        my $R2 = $f[-1];
        my $output= "$snp_1\t$snp_2\t$R2";
        foreach my $k (sort keys %hash1 ){
            my $bind = $hash1{$k};
            my @t= split/\t/,$k;
            my $chr =$t[0];
            my $start =$t[1];
            my $end = $t[2];
            if ($chr1 == $chr){
                if ($snp1>= $start  && $snp1 <=$end && $snp2 >= $start && $snp2 <= $end ){
                    print $O1 "$output\t$bind\n";
                }
                else{
                    unless(exists $hash2{$output}){
                        $hash2{$output} =1;
                        print $O2 "$output\n";
                    }
                }
            }
        }
    }
}
