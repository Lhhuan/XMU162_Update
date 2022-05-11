#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

use List::MoreUtils ':all';
my $f4 = "../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz";
# my $f4 = "../output/01_hotspot_target_gene_reactomeFI.bed.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

while(<$I4>)
{
    chomp;
    my @f= split/\t/;
    unless(/^hotspot_chr/){
        my $egene= $f[3];
        my $co=$f[-1];
        my @t=split/;/,$co;
        my $num=@t;
        print "$egene\t$num\n";
    }
}