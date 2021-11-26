#将"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz" 中的../output/01_top_QTL_eQTL_mQTL_miQTL.txt.gz 提取出来，
#得../output/02_hotspot_in_1kg_phase3_v5_eur.vcf.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my $f1 = "../output/01_top_QTL_eQTL_mQTL_miQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/02_hotspot_in_1kg_phase3_v5_eur.vcf.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/02_hotspot_out_1kg_phase3_v5_eur.vcf.gz";
open my $O2, "| gzip >$fo2" or die $!;


my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/SNP/){
        my @f = split/\t/;
        my $SNP =$f[0];
        $hash1{$SNP}=1;
    }
}



while(<$I2>)
{
    chomp;
    if(/#/){
        print $O1 "$_\n";
    }
    else{
        my @f = split/\s+/;
        my $ID =$f[2];
        if (exists $hash1{$ID}){
            print $O1 "$_\n";
            $hash2{$ID}=1;
        }
    }
}

foreach my $k(sort keys %hash1){
    unless(exists $hash2{$k}){
        print $O2 "$k\n";
    }
}