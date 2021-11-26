#将/home/chaoqun/phase3/integrated_call_samples_v3.20130502.ALL.panel 中的EUR提出来，得EUR_sample_list.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "new_file.vcf";
# my $f1 = "head_chr9";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# my $fo1 = "EUR_sample_list.txt";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


while(<$I1>)
{
    chomp;
    if (/^#CHROM/){
        my @f =split/\s+/;
        my $number =@f;
        print "$number\n";
    }
}
