#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my @a=(2,8,5,6);
my @sorted_vs = sort {$a <=> $b} @a;
# print 

my $cc  =join("\t",@sorted_vs);
print "$cc\n";