#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
my $f1 = "03_entrezgene_symbol_ensembl_pos.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "01_G16808_S85825_get_gene_file_list.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my $fo2 = "04_entrezgene_ensembl.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "entrezgene\tensembl\tChr\tstart\tend\n";

my %hash1;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^entrezgene/){
       my $entrezgene =$f[0];
       my $ensembl =$f[1];
       my $chr=$f[-3];
       my $start =$f[-2];
       my $end = $f[-1];
       my $v ="$ensembl\t$chr";
       unless($ensembl =~/NULL/){
           $hash1{$entrezgene}=$v;    
           print $O2 "$entrezgene\t$ensembl\t$chr\t$start\t$end\n";
       }

    }
}

my @files=<$I2>;
foreach my $file_entrezgene(@files){
    chomp($file_entrezgene);
    if(exists $hash1{$file_entrezgene}){
        my $v = $hash1{$file_entrezgene};
        my @t=split/\t/,$v;
        my $file_ensembl =$t[0];
        my $file_chr = $t[1];
        # print "$v\n";
        my $f3 = "./G16808_S85825/${file_entrezgene}";
        open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
        my $fo1 = "./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}.txt.gz"; #
        # open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
        open my $O1, "| gzip >$fo1" or die $!;
        my $fo2 = "./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}_same_pos.txt.gz"; #
        open my $O2, "| gzip >$fo2" or die $!;
        while(<$I3>)
        {
            chomp;
            my @f= split/\t/;
            my $entrezgene =$f[0];
            # print "$.\n";
            if($. <=300 && exists $hash1{$entrezgene}){
                # print "$entrezgene\n";
                my $vs = $hash1{$entrezgene};
                my @tt = split/\t/,$vs;
                my $ensg =$tt[0];
                my $ensg_chr= $tt[1];
                print $O1 "$ensg\t$entrezgene\n";
                if($file_chr eq $ensg_chr){
                    print $O2 "$ensg\t$entrezgene\n";
                }
            }
        }
        print "$file_entrezgene\n";
    }
}
