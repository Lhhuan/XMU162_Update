#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $dir = "./G16808_S85825";
opendir (DIR, $dir) or die "can't open the directory!";
my @dir = readdir DIR;
my $fo1 = "01_G16808_S85825_get_gene_file_list.txt.gz"; #
open my $O1, "| gzip >$fo1" or die $!;

foreach my $file(@dir){
    if($file=~/^\d+$/){
        print $O1 "$file\n";
    }
}