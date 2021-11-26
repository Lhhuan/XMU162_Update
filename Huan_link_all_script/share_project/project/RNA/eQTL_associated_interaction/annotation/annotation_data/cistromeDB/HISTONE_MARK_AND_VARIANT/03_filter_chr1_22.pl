#过滤出merge_pos_info_narrow_peak_sorted.bed.gz中的chr1_22,得03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

my $f1 = "merge_pos_info_narrow_peak_sorted.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;
my %hash2;

for (my $i=1;$i<23;$i++){
    my $k = "chr${i}";
    $hash1{$k}=1;
}


while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];
    if (exists $hash1{$chr}){
        print $O1 "$_\n";
    }
}
