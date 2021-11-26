#将../data/mat_rank_gs_split_sig.txt中没有转成 symbol用原来的gene表示，
#得../07_final_all_gene_from_jintao_symbol.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 ="../output/06_all_gene_from_jintao_symbol.txt";
my $f2 ="../data/mat_rank_gs_split_sig.txt";
my $fo1 ="../output/07_final_all_gene_from_jintao_symbol.txt"; 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


my(%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7);

print $O1 "Gene\tSymbol\tKinase\n";
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
    unless(/^gene/){
        my @f= split /\t/;
        my $gene  = $f[0];
        my $kinase = $f[-1];
        if (exists $hash1{$gene}){
            my $symbol =$hash1{$gene};
            my $output = "$gene\t$symbol\t$kinase";
            unless(exists $hash2{$output}){
                $hash2{$output}=1;
                print $O1 "$output\n";
            }
            
        }
        else{  #没能转换成symbol的gene 用其原来的gene作为symbol，
            my $gene1 = $gene;
            $gene1 =~ s/\s+/_/g; #为了后面排序方便，将作为symbol的gene名字中的空格替换为_
            my $output = "$gene\t$gene1\t$kinase";
            unless(exists $hash2{$output}){
                $hash2{$output}=1;
                print $O1 "$output\n";
            }
        }
    }
}
