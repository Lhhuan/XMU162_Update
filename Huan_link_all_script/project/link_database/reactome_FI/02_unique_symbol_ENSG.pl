#unique ./output/01_symbol_ENSG.txt,得"./output/02_query_gene_entrezgene_symbol_ensembl.txt" 和"./output/query_gene_ensembl.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "./output/01_symbol_ENSG.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "./output/02_query_gene_entrezgene_symbol_ensembl.txt"; #
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "./output/query_gene_ensembl.txt"; #
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";

print $O1 "query_gene\tentrezgene\tsymbol\tensembl\n";
print $O2 "query_gene\tensembl\n";




my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^query/){
        my $gene =$f[0];
        my $ensg = $f[-3];
        my $Entrez = $f[-2];
        my $SYMBOL =$f[-1];
        $ensg =~ s/\s+//g;
        $ensg =~ s/list//g;
        $ensg =~ s/gene=//g;
        $ensg =~ s/,.*//g;
        $ensg =~ s/\(//g;
        $ensg =~ s/\)//g;
        $ensg =~ s/c//g;
        $ensg =~ s/"//g;
        my $v ="$Entrez\t$SYMBOL\t$ensg";
        push @{$hash1{$gene}},$v;
        push @{$hash2{$gene}},$ensg;
    }
}


foreach my $k(sort keys %hash1){
    unless($k=~/ENSG/){
        my @v1s =@{$hash1{$k}};
        my $v1 = $v1s[0];
        print $O1 "$k\t$v1\n";
        my @v2s=@{$hash2{$k}};
        my $v2=$v2s[0];
        print $O2 "$k\t$v2\n";
    }

}
