#对合并enhancer annotation"../output/nodes_annotation/EnhancerAtlas_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"，"../output/nodes_annotation/JEME_enhancer.bed.gz"信息得"../output/nodes_annotation/Adjust_${marker}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my $marker = "Enhancer";
my $f1 = "../output/nodes_annotation/EnhancerAtlas_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/nodes_annotation/JEME_enhancer.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n");
my $fo1 = "../output/nodes_annotation/Adjust_${marker}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;

# print $O1 "hchr_hstart_hend\tehchr_ehstart_ehend:cell_line:Source\n";
print $O1 "Hotspot\tEnhancer\n";
my %hash1;
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $h_chr=$f[0];
    my $h_start =$f[1];
    my $h_end =$f[2];
    my $f_pos =join("_",@f[3..5]);
    my $cell_line =$f[6];
    my $k= join("_",@f[0..2]);
    my $v="$f_pos:$cell_line:EnhancerAtlas";
    push @{$hash1{$k}},$v;
}

while(<$I2>)
{
    chomp;
    my @f =split/\t/;
    my $h_chr=$f[0];
    my $h_start =$f[1];
    my $h_end =$f[2];
    my $f_pos =join("_",@f[4..6]);
    my $cell_line =$f[-2];
    my $source = $f[-1];
    my $k= join("_",@f[0..2]);
    my $v="$f_pos:$cell_line:$source";
    push @{$hash1{$k}},$v;
}
foreach my $k(sort keys %hash1){
    my @vs= @{$hash1{$k}};
    @vs = uniq(@vs);
    my $v=join(";",@vs);
    print $O1 "$k\t$v\n";
}
close($O1);
        