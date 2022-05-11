# 调整../output/nodes_annotation/hotspot_tissue_anno.bed.gz 的格式得../output/nodes_annotation/Adjust_hotspot_tissue_anno.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';


my $f1 = "../output/nodes_annotation/hotspot_tissue_anno.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/nodes_annotation/Adjust_hotspot_tissue_anno.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
    
my %hash1;

print $O1 "Hotspot\teQTL_Source\n"; #postion:cell_line
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $h_chr=$f[0];
    my $h_start =$f[1];
    my $h_end =$f[2];
    my $f_pos =join("_",@f[3..5]);
    my $tissue =$f[7];
    my $k= join("_",@f[0..2]);
    my $v = "$f_pos:$tissue";
    push @{$hash1{$k}},$v;
}
foreach my $k(sort keys %hash1){
    my @vs= @{$hash1{$k}};
    @vs = uniq(@vs);
    my $v=join(";",@vs);
    print $O1 "$k\t$v\n";
}
        