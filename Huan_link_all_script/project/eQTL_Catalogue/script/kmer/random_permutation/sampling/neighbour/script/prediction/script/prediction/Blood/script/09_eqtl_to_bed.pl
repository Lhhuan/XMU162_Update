#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/GTEx/GTEx_ge_blood.all.tsv.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo2 = "../output/09_eqtl.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# print $O2 "SNP_chr\tSNP_pos\tPvalue\tegene\ttissue\tqtl_group\n";



while(<$I2>)
{
    chomp;
    unless(/^molecular_trait_id/){
        my @ss = split/\t/;
        my $chr = $ss[1];
        my $pos= $ss[2];
        my $ref= $ss[3];
        my $alt= $ss[4];
        # my $variant =$ss[5];
        my $pvalue =$ss[8];
        my $chrpos ="$chr\t$pos";
        my $variant = "$chrpos\t$ref\t$alt";
        my $gene_id = $ss[-3];
        # print "$variant\n$chrpos\n";
        my $start = $pos-1;
        my $end=$pos+1-1;
        print $O2 "chr${chr}\t$start\t$end\t$pvalue\t$gene_id\n";
    }
}
close($O2);
system "zless $fo2 |sort -k1,1 -k2,2n |gzip >../output/09_eqtl_sorted.bed.gz";