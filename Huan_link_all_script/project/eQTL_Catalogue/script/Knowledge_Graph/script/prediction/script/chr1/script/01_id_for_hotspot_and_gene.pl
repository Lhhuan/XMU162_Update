#利用"../../../../../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz" 为chr1 hotspot和gene 建立index分别得"../output/01_hotspot_idx.txt.gz","../output/01_gene_idx.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';


my $f1 = "../../../../../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz";
# my $f1 = "1234.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/01_hotspot_idx.txt.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Hotspot\thotspot_idx\n";
my $fo2 = "../output/01_gene_idx.txt.gz"; #
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "Gene\tGene_idx\n";
my (%hash1,%hash2);



while(<$I1>)
{
    chomp;
    unless(/^hotspot_chr/){
        my @f =split/\t/;
        my $hotspot_chr = $f[0];
        if($hotspot_chr =~/^chr1$/){
            my $hotspot=join("_",@f[0..2]);
            $hash1{$hotspot}=1;
            my $egene =$f[3];
            $hash2{$egene}=1;
            my $ReactomeFI_gene =$f[6];
            my $egene_co_expression_gene =$f[7];
            unless($ReactomeFI_gene =~/NA/){
                my @ts= split/;/,$ReactomeFI_gene;
                foreach my $t(@ts){
                    $hash2{$t}=1;
                }
            }
            unless($egene_co_expression_gene =~/NA/){
                my @ts= split/;/,$egene_co_expression_gene;
                foreach my $t(@ts){
                    $hash2{$t}=1;
                }            
            }
        }
    }
}

my @keys1 =keys %hash1;

foreach my $i(0..$#keys1){
    print $O1 "$keys1[$i]\t$i\n";
}

my @keys2 =keys %hash2;
print "$#keys2\n";
my $start=$#keys1+1;
my $end=$start+$#keys2;
foreach my $j($start..$end){
    my $true=$j-$#keys1-1;
    print $O2 "$keys2[$true]\t$j\n";
}