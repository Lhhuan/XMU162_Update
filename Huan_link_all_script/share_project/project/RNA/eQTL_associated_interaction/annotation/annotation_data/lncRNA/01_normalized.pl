#将NONCODEv5_hg19.lncAndGene.bed.gz中需要的列提出来，得NONCODEv5_hg19.lncAndGene_need.bed.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "NONCODEv5_hg19.lncAndGene.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "NONCODEv5_hg19.lncAndGene_need.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "chr\tstart\tend\n";


while(<$I1>)
{
    chomp;
    my @f = split/\s+/;
    my $chr = $f[0];
    my $start =$f[1];
    my $end =$f[2];
    print $O1 "$chr\t$start\t$end\n";
}
