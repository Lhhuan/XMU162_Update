#将01_hsa_hg19.gff3 转换为01_hsa_hg19.bed.gz，bed为0-based文件

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $f1 = "01_hsa_hg19.gff3";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "01_hsa_hg19.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "chr\tstart\tend\tmiRNA_type\tStrand\tinfo\n";

while(<$I1>)
{
    chomp;
    unless(/^Chr/){
        my @f = split/\s+/;
        my $chr = $f[0];
        my $miRNA_type =$f[2];
        my $start =$f[3];
        my $end =$f[4];
        my $Strand =$f[6];
        my $info = $f[-1];
        my $new_start = $start -1;
        my $output = "$chr\t$new_start\t$end\t$miRNA_type\t$Strand\t$info";
        print $O1 "$output\n";
    }
}