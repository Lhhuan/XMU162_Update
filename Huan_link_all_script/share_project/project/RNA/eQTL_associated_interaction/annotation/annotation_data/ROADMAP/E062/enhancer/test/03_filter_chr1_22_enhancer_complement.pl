#将01_enhancer_complement.bed中的chr1 和chr22提取出来，得01_enhancer_complement_chr1_22.bed
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;

# my $f1 = "../output/Whole_Blood/Cis_eQTL/NHP/interval_${j}_chr1.txt.gz";
my $f1 = "01_enhancer_complement.bed";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo3 = "01_enhancer_complement_chr1_22.bed";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
# open my $O3, "| gzip >$fo3" or die $!;
my %hash1;
while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];

    if ($chr=~/\bchr1\b/ || $chr=~/\bchr22\b/){
        print $O3 "$_\n";
    }
}
