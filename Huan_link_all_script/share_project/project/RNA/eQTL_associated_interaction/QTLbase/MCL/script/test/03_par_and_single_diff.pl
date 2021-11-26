#将用 ../data/cytoBand.txt.gz对compute-0-6下的/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz 
#进行annotation,得../output/03_annotation_ld_band.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my $f1 = "/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld_x000_validate.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/state/partition1/huan/tmp_output/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld_x000_annotation_bind.txt.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "/state/partition1/huan/03_single_in_par.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $header ="SNP_A\tSNP_B\tR2";
print $O1 "$header\tbind\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/^SNP_A/){
        $hash1{$_}=1;
    }
}



while(<$I2>)
{
    chomp;
    unless (/SNP_A/){
       my $f = $_;
       if (exists $hash1{$f}){
           print $O1 "$_\n";
       }
       else{
           print "$_\n";
       }
    }
}
