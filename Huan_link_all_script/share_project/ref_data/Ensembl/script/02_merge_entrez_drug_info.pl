#将../output/Hg19_reference_geneome_data_from_ensembl_BioMart.txt 和../output/01_transfrom_ensg_to_entrez.txt merge 到一起，
#得../output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/01_transfrom_ensg_to_entrez.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../output/Hg19_reference_geneome_data_from_ensembl_BioMart.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 = "../output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";

my %hash1;

while(<$I1>) #去掉重复的header
{
    chomp;
    unless(/^query/){
        my @f =split/\t/;
        my $ensg =$f[0];
        my $entrezgene =$f[1];
        my $symbol =$f[2];
        my $v = "$entrezgene\t$symbol";
        $hash1{$ensg}=$v;
    }
}

while(<$I2>) #去掉重复的header
{
    chomp;
    if(/^Gene/){
        print $O1 "$_\tEntrezgene\tSymbol\n";
    }
    else{
        my @f =split/\t/;
        my $ensg =$f[0];
        if (exists $hash1{$ensg}){
            my $entrezgene =$hash1{$ensg};
            print $O1 "$_\t$entrezgene\n";
        }
        else{
            print $O1 "$_\tNA\tNA\n";
        }
    }
}
    