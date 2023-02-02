#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
my $fo1 = "../../../output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../../../output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz";
my $fo3 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_gwas.bed.gz";
my $fo4 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_gwas_metaid.bed.gz";
open my $O4, "| gzip >$fo4" or die $!;
print $O4 "CHR\tstart\tend\tcluster\tGWAS_chr\tGWAS_start\tGWAS_end\tmeta_id\n";

my %hash1;
my $f1 = "../../../output/figure/whole_genome/11_whole_genome_leiden_pca3_k50_resolution2e-05.txt";
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
# system "bedtools intersect -a $fo2 -b /share/data0/QTLbase/huan/gwas_catalogy/from_guojintao_gwas_catalogue_sort.bed.gz -wa -wb |gzip >$fo3 ";

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
    my $gwas_chr= $f[4];
    my $gwas_start=$f[5];
    my $gwas_end=$f[6];
    my $info=$f[-1];
    my @t=split/\;/,$info;
    my $meta_id=$t[-3];
    my $output= "$chr\t$start\t$end\t$cluster\t$gwas_chr\t$gwas_start\t$gwas_end\t$meta_id";
    unless(exists $hash1{$output}){
        $hash1{$output}=1;
        print $O4 "$output\n";
    }
}

my $f3 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GWAS/random_hotspot_from_whole_genome_causaldb.bed.gz";
system "bedtools intersect -a /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/random_hotspot_from_whole_genome.bed.gz -b /share/data0/QTLbase/huan/CAUSALdb/output/02_EUR_causaldb_variant_hg38_sorted.bed.gz -wa -wb |gzip > $f3 ";


my $fo5 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GWAS/17_whole_genome_random_causaldb_gwas_metaid.bed.gz";
open my $O5, "| gzip >$fo5" or die $!;
print $O5 "CHR\tstart\tend\tGWAS_chr\tGWAS_start\tGWAS_end\tmeta_id\n";


open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); 
my %hash2;
while(<$I3>)
{
    chomp;
    my @f= split/\t/;
    my $chr = $f[0];
    my $start  =$f[1];
    my $end = $f[2];
    my $gwas_chr= $f[3];
    my $gwas_start=$f[4];
    my $gwas_end=$f[5];
    my $info=$f[-1];
    my @t=split/\;/,$info;
    my $meta_id=$t[-3];
    my $output= "$chr\t$start\t$end\t$gwas_chr\t$gwas_start\t$gwas_end\t$meta_id";
    unless(exists $hash2{$output}){
        $hash2{$output}=1;
        print $O5 "$output\n";
    }
}