#"/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 注释ensg位置得../output/09_hotspot_0.05egene_pos.bed.gz
#与TAD /home/huanhuan/project/link_database/ENCOED/output/hg38/01_merge_TAD_sample_sorted.bed.gz intersect得"../output/edges_annotation/09_TAD_hotspot.bed.gz"， 将gene位置放前面得"../output/edges_annotation/09_Adjust_TAD_hotspot.bed.gz"，用gene 位置与TAD intersect得"../output/edges_annotation/09_TAD_egene.bed.gz"， 判断hotspot-egene在相同TAD得"../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my (%hash1,%hash2);

my $fo1 = "../output/09_hotspot_0.05egene_pos.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Hotspot_gene\t${type}\n";
my $f1= "/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2= "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
while(<$I1>)
{
    chomp;
    unless(/Gene/){
        my @f=split/\t/;
        my $ensg =$f[0];
        my $start =$f[1];
        my $end=$f[2];
        my $chr=$f[3];
        unless($chr=~/CHR|GL|KI|MT|X|Y/){
            my $v= "chr${chr}\t$start\t$end";
            $hash1{$ensg}=$v;
        }
    }
}

while(<$I2>)
{
    chomp;
    my @f=split/\t/;
    my $ensg=$f[3];
    if(exists $hash1{$ensg}){
        my $gene_pos = $hash1{$ensg};
        print $O1 "$_\t$gene_pos\n";
    }
}

close($O1);
my $f3 = "../output/edges_annotation/09_TAD_hotspot.bed.gz";
system "bedtools intersect -a $fo1 -b /home/huanhuan/project/link_database/ENCOED/output/hg38/01_merge_TAD_sample_sorted.bed.gz -wa -wb |gzip >$f3";


open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
my $fo2 = "../output/edges_annotation/09_Adjust_TAD_hotspot.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
while(<$I3>)
{
    chomp;
    my @f=split/\t/;
    my $hotspot_p =join("\t",@f[0..2]);
    my $egene =$f[3];
    my $gene_p =join("\t",@f[4..6]);
    my $TAD_p = join("\t",@f[7..9]);
    my $tissue =$f[-1];
    print $O2 "$gene_p\t$hotspot_p\t$egene\t$TAD_p\t$tissue\n";
    # my $k = 
}
close($O2);
my $f4 = "../output/edges_annotation/09_TAD_egene.bed.gz";
system "bedtools intersect -a $fo2 -b /home/huanhuan/project/link_database/ENCOED/output/hg38/01_merge_TAD_sample_sorted.bed.gz -wa -wb |gzip >$f4";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

my $fo3 = "../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;
while(<$I4>)
{
    chomp;
    my @f=split/\t/;
    my $gene_p =join("\t",@f[0..2]);
    my $hotspot_p =join("\t",@f[3..5]);
    my $egene =$f[6];
    my $TAD_p1 = join("\t",@f[7..9]);
    my $tissue1 =$f[10];
    my $TAD_p2 = join("\t",@f[11..13]);
    my $tissue2 =$f[14];
    my $k1= "$TAD_p1\t$tissue1";
    my $k2= "$TAD_p2\t$tissue2";
    if($k1 eq $k2){
    print $O3 "$gene_p\t$hotspot_p\t$egene\t$TAD_p1\t$tissue1\n";
    }
}