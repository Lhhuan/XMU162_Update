#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $region_start = 144275000;
my $region_end = $region_start+ 10e+6;
print "chr1\t$region_start\t$region_end\n";