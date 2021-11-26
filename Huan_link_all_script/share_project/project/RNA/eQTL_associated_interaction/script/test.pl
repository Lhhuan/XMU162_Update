#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

# my $gene = "BSCL2///HNRPUL2";
# my $ddd= "///";
# $gene =~ s/$ddd/,/g;
# print "$gene\n";

my $line1s =readpipe("wc -l /home/huanhuan/project/RNA/eQTL_associated_interaction/script/unique_gene_name.txt");
my @c = split/\s+/,$line1s;
my $line = $c[0];
print "$line\n";