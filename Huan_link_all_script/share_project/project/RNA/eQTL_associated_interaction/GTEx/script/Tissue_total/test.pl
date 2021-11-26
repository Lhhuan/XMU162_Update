#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @ds = (1,4,5,5e-8,1e-5,10E-7);
foreach my $d(@ds){
    if($d <= 1E-5){
        print "$d\n";
    }
}


# my $a = 1e-7;
# my $b = 1e-5;
# my $c = $b-$a;
# print "$c\n";
