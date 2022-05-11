#调整"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"格式得 "../output/01_all_hotspot.csv.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/01_all_hotspot.csv.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
# my $fo2 = "../output/train/01_hotspot.csv.gz"; #
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# open my $O2, "| gzip >$fo2" or die $!;
print $O1 "Source node name,Target node name,Source node type,Target node type,Relationship type\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $hotspot = join("_",@f[0..2]);
    my $egene =$f[-1];
    print $O1 "$egene,$hotspot,gene,hotspot,risk_gene\n";
}




