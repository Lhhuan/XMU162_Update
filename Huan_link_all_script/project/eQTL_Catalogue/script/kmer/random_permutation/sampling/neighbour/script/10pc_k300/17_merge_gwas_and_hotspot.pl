#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
my $fo1 = "../output/figure/09_chr/17_chr1_louvain_pca10_k300_hotspot.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/figure/09_chr/17_chr1_louvain_pca10_k300_hotspot_sorted.bed.gz";
my $fo3 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/GWAS/17_chr1_louvain_pca10_k300_hotspot_gwas.bed.gz";
my $fo4 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/GWAS/17_chr1_louvain_pca10_k300_hotspot_gwas_metaid.bed.gz";
open my $O4, "| gzip >$fo4" or die $!;
print $O4 "CHR\tstart\tend\tcluster\tmeta_id\n";

my %hash1;
my $f1 = "../output/figure/09_chr/11_chr1_louvain_pca10_k300.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

while(<$I1>)
{
    chomp;
    unless(/^PC_1/){
        my @f= split/\t/;
        my $cluster=$f[-2];
        my $hotspot=$f[-1];
        my @t=split/\:/,$hotspot;
        my $chr = $t[0];
        my @ss=split/-/,$t[1];
        my $start = $ss[0];
        my $end =$ss[1];
        print $O1 "$chr\t$start\t$end\t$cluster\n";
    }
}

close($O1);
system"zless $fo1 |sort -k1,1 -k2,2n |gzip >$fo2";
system "bedtools intersect -a $fo2 -b /share/data0/QTLbase/huan/CAUSALdb/output/02_EUR_causaldb_variant_hg38_sorted.bed.gz -wa -wb |gzip >$fo3 ";

my $f2 = $fo3;
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); 

while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    my $chr = $f[0];
    my $start  =$f[1];
    my $end = $f[2];
    my $cluster=$f[3];
    my $info=$f[-1];
    my @t=split/\;/,$info;
    my $meta_id=$t[-3];
    my $output= "$chr\t$start\t$end\t$cluster\t$meta_id";
    unless(exists $hash1{$output}){
        $hash1{$output}=1;
        print $O4 "$output\n";
    }
}