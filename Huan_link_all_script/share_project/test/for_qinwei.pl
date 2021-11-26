#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "/home/qinwei/tcgabarcode_AML.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");

my $f2 = "/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/mc3.v0.2.8.PUBLIC.maf.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); 

my $fo1 = "filter_mc3.v0.2.8.PUBLIC.maf.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;


while(<$I1>)
{
    chomp;
    my $a = $_;
    $a =~ s/\s+//g;
    $hash1{$a}=1;

}

while(<$I2>)
{
    chomp;
    my @f = split/\t/;
    my $sample = $f[15];
    if(/^Hugo_Symbol/){
        print $O1 "$_\n";
        # print "$sample\n";
    }
    else{
        my @t = split/\-/,$sample;
        my $new_sample_format = join("-",@t[0..2]);
        $new_sample_format =~ s/\s+//g;
        # print "$sample\n";
        # print "$new_sample_format\n";
        if (exists $hash1{$new_sample_format}){
            print $O1 "$_\n";
        }
    }
}