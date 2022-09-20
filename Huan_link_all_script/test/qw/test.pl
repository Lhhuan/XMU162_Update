#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "sample121_freq.frq";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $fo1 = "output.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


while(<$I1>)
{
    chomp;
    my @f=split/\t/;
    for (my $i=0;$i<6;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    unless(defined $f[$i]){
    $f[$i] = "NA";
    }
    unless($f[$i]=~/\w/){$f[$i]="NA"} #对文件进行处理，把所有定义的没有字符的都替换成NULL
    }
    my $t =join("\t",@f[0..5]);
    print $O1 "$t\n";
}


