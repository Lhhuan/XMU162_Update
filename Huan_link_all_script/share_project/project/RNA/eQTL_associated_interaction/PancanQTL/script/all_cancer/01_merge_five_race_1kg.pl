#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $f2 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my $f3 = "/share/data0/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz";
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

my $f4 = "/share/data0/1kg_phase3_v5_hg19/SAS/1kg.phase3.v5.shapeit2.sas.hg19.all.SNPs.vcf.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

my $f5 = "/share/data0/1kg_phase3_v5_hg19/AMR/1kg.phase3.v5.shapeit2.amr.hg19.all.SNPs.vcf.gz";
open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件

my $f6 = "/share/data0/1kg_phase3_v5_hg19/AFR/1kg.phase3.v5.shapeit2.afr.hg19.all.SNPs.vcf.gz";
open( my $I6 ,"gzip -dc $f6|") or die ("can not open input file '$f6' \n"); #读压缩文件
#---------------------------mkdir
#---------------
my $fo1 = "../../output/five_race_1kg_snp.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "SNP_chr\tSNP_pos\n";


my %hash1;
while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\n";
            $hash1{$k}=1;
        }
    }
}
print "111\n";
while(<$I3>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\n";
            $hash1{$k}=1;
        }
    }
}

print "2\n";
while(<$I4>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\n";
            $hash1{$k}=1;
        }
    }
}

print "3\n";
while(<$I5>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\n";
            $hash1{$k}=1;
        }
    }
}

print "4\n";
while(<$I6>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\n";
        }
    }
}
print "5\n";