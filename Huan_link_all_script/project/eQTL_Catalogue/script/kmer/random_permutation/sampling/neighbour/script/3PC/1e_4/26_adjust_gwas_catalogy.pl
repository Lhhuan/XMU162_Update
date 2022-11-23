#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;
my $f1 = "/share/data0/QTLbase/huan/gwas_catalogy/gwasCatalog.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 

my $fo2 = "/share/data0/QTLbase/huan/hic/4DNFIFLJLIS5_5000.ginteractions_cross_chr_adjust.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $bin=$f[0];
    my $chr=$f[1];
    my $start =$f[2];
    my $end=$f[3];
    my $rsid=$f[4];
    my $trait =$f[10];
    my $output="$chr\t$start\t$end\t$rsid\t$trait\n";
    print $O2 "$output\n";
}
close($O2);
my $sorted_fo2 ="/share/data0/QTLbase/huan/gwas_catalogy/gwasCatalog_adjust_sorted.bed.gz";
system "zless $fo2 |sort -k1,1 -k2,2n |gzip >$sorted_fo2";

system "bedtools intersect -a ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -b $sorted_fo2 -wa -wb |gzip >../../../output/figure/whole_genome/3pca_1e_4/GWAS/26_gwas_catalogy_cluster.bed.gz";







