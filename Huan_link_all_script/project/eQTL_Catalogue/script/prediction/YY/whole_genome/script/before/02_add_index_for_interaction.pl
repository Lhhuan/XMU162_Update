#用"../output/01_hotspot_idx.txt.gz" 和"../output/01_gene_idx.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"添加idx,得"../output/02_hotspot_egene_idx.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

my $f1 = "../output/01_hotspot_idx.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/01_gene_idx.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $f3 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz";
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件


my $fo1 = "../output/02_hotspot_egene_idx.txt.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Hotspot_idx\tEgene_idx\n";

my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    unless(/^Hotspot/){
        my @f =split/\t/;
        my $hotspot  =$f[0];
        my $idx=$f[1];
        $hash1{$hotspot}=$idx;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^Gene/){
        my @f =split/\t/;
        my $Gene  =$f[0];
        my $idx=$f[1];
        $hash1{$Gene}=$idx;
    }
}

while(<$I3>)
{
    chomp;
    unless(/^hotspot_chr/){
        my @f =split/\t/;
        if($f[0]=~/^chr1$/){
            my $hotspot=join("_",@f[0..2]);
            my $egene =$f[3];
            if(exists $hash1{$hotspot}){
                my $hotspot_idx = $hash1{$hotspot};
                if(exists $hash1{$egene}){
                    my $egene_idx =$hash1{$egene};
                    print $O1 "$hotspot_idx\t$egene_idx\n";
                }
            }
        }
    }
}
