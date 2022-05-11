#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
my $f1 = "03_entrezgene_symbol_ensembl.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "01_G16808_S85825_get_gene_file_list.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my $fo2 = "03_entrezgene_ensembl.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "entrezgene\tensembl\n";

my %hash1;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^entrezgene/){
       my $entrezgene =$f[0];
       my $ensembl =$f[1];
       unless($ensembl =~/NULL/){
           $hash1{$entrezgene}=$ensembl;    
           print $O2 "$entrezgene\t$ensembl\n";
       }

    }
}

my @files=<$I2>;
foreach my $file_entrezgene(@files){
    chomp($file_entrezgene);
    if(exists $hash1{$file_entrezgene}){
        my $file_ensembl = $hash1{$file_entrezgene};
        my $f3 = "./G16808_S85825/${file_entrezgene}";
        open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
        my $fo1 = "./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}.txt.gz"; #
        # open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
        open my $O1, "| gzip >$fo1" or die $!;
        while(<$I3>)
        {
            chomp;
            my @f= split/\t/;
            my $entrezgene =$f[0];
            # print "$.\n";
            if($. <=300 && exists $hash1{$entrezgene}){
                # print "$entrezgene\n";
                my $ensg = $hash1{$entrezgene};
                print $O1 "$ensg\t$entrezgene\n";
                # print  "$ensg\t$entrezgene\n";
            }
        }
        print "$file_entrezgene\n";
    }
}
