#unique "02_G16808_S85825_entrezgene_symbol_ensembl.gene.txt" 的ENSG,得03_entrezgene_symbol_ensembl.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "02_G16808_S85825_entrezgene_symbol_ensembl.gene.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "03_entrezgene_symbol_ensembl.txt.gz"; #
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "entrezgene\tensembl\tsymbol\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^query/){
        my $Entrez =$f[0];
        my $ensg = $f[-3];
        my $SYMBOL =$f[-2];
        $ensg =~ s/\s+//g;
        $ensg =~ s/list//g;
        $ensg =~ s/gene=//g;
        $ensg =~ s/,.*//g;
        $ensg =~ s/\(//g;
        $ensg =~ s/\)//g;
        $ensg =~ s/c//g;
        $ensg =~ s/"//g;
        # my $v ="$SYMBOL\t$ensg";
        print $O1 "$Entrez\t$ensg\t$SYMBOL\n";
    }
}