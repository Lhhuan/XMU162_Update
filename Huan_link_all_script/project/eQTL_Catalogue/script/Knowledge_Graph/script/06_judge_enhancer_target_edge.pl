#筛选 ../output/edges_annotation/enhancer_target_anno.bed.gz中enhancer-target gene和 egene是否是相同基因，得../output/edges_annotation/success_enhancer_target_anno.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';


my $f1 = "../output/edges_annotation/enhancer_target_anno.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/edges_annotation/success_enhancer_target_anno.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
    
my %hash1;

print $O1 "Hotspot_egene\tenhancer_target_source\n"; #postion:cell_line
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $h_chr=$f[0];
    my $h_start =$f[1];
    my $h_end =$f[2];
    my $egene = $f[3];

    my $f_pos =join("_",@f[4..6]);
    my $enhancer_target_gene= $f[-4];
    my $tissue =$f[-2];
    if($egene eq $enhancer_target_gene){
        # print "$_\n";
        my $k="${h_chr}_${h_start}_${h_end}:${egene}";
        my $v = "$f_pos:$tissue";
        push @{$hash1{$k}},$v;
    }
}
foreach my $k(sort keys %hash1){
    my @vs= @{$hash1{$k}};
    @vs = uniq(@vs);
    my $v=join(";",@vs);
    print $O1 "$k\t$v\n";
}
        