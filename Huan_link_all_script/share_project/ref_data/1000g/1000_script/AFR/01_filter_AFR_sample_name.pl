#将/home/chaoqun/phase3/integrated_call_samples_v3.20130502.ALL.panel 中的AFR提出来，得AFR_sample_list.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/chaoqun/phase3/integrated_call_samples_v3.20130502.ALL.panel";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "AFR_sample_list.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


while(<$I1>)
{
    chomp;
    unless(/^sample/){
        my @f =split/\t/;
        my $sample = $f[0];
        my $super_pop =$f[2];
        if ($super_pop =~/AFR/){
            print $O1 "$sample\n";
        }
    }
}
