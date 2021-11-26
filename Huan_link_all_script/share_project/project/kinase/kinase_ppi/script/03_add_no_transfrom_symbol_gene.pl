#../output/all_network_gene.txt 中 没能转换成symbol的gene 用其原来的gene代替，
#得../output/03_final_all_network_gene_symbol.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 ="../output/02_all_network_gene_symbol.txt";
my $f2 ="../output/all_network_gene.txt";
my $fo1 ="../output/03_final_all_network_gene_symbol.txt"; 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


my(%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7);

print $O1 "Gene\tSymbol\n";
while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^query/){
        my $gene = $f[0];
        my $symbol = $f[1];
        $hash1{$gene}=$symbol;
    }
}

while(<$I2>)
{
    chomp;
    my @f= split /\t/;
    my $gene  = $f[0];
    if (exists $hash1{$gene}){
        my $symbol =$hash1{$gene};
        print $O1 "$gene\t$symbol\n";
    }
    else{  #没能转换成symbol的gene 用其原来的gene作为symbol，
        my $gene1 = $gene;
        $gene1 =~ s/\s+/_/g; #为了后面排序方便，将作为symbol的gene名字中的空格替换为_
        print $O1 "$gene\t$gene1\n";
    }
}
