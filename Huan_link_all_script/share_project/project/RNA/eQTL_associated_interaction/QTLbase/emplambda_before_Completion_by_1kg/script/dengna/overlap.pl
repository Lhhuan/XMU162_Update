#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中按照位置，两两取交集，得../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
# use List::Util qw/max min/;

my %hash2;
my $f1 = "gene1.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "gene2.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 = "overlap_gene.txt";
# open my $O1, "| gzip >$fo1" or die $!;
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# print $O1 "${QTL1}_emplambda\t${QTL2}_emplambda\tpos\tchr\n";
# print "$cp_format1\n";

while(<$I1>)
{
    chomp;  
    my @f = split/\t/;
    my $gene1 =$f[0]; 
    $gene1 =~ s/^s+//g;
    $gene1 =~ s/s+&//g;
    $hash2{$gene1}=1;
    if($gene1 =~/RRP7A/){
        print "$_\t222\n";
    }
    

}

while(<$I2>)
{
    chomp;  
    my @f = split/\t/;
    my $gene2 =$f[0]; 
    $gene2 =~ s/^s+//g;
    $gene2 =~ s/s+&//g;
    # print "$gene2\n";
    if(exists $hash2{$gene2}){
        print $O1 "$gene2\n";
        # print "$gene2\n";
    }
    if($gene2 =~/RRP7A/){
        print "$_\n";
    }

}