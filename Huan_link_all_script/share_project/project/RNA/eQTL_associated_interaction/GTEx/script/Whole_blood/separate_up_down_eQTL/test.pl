#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my @files = (1,2,4,5,6,3);
my @sorted_files = sort {$b <=> $a} @files;
print "@sorted_files\n";