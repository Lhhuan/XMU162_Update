#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


# my $a =int(23.5);
# my $b = 24 % 2;
my @cc=(1,23,5,6,7);
my @cc1 = sort { $a <=> $b } @cc;
# my @cc1 = sort { $b <=> $a } @cc;
print "@cc1\n";


# my @a1=(4,5,1,3,6,2,10);
# my $ordered = join ",",sort {$a<=>$b} @a1;
# print "ordered $ordered\n"