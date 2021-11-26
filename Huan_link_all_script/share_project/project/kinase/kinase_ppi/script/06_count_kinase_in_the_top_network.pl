#根据../output/01_transfrom_kinsae_gene_to_symbol.txt和../output/RWR/rwr_result_top/ 计算kinase in the top, 
#and kinase out of top, isn't no-kinase in the top and isn't no-kinase out of the top,得文件../output/06_fisher_need_data.txt
#得不在网络中的kinase文件../output/06_kinase_not_in_the_network.txt
#得在网络中的kinase文件 ../output/06_kinase_in_the_network.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use File::Copy;

my $f1 ="../output/03_final_all_network_gene_symbol.txt";
my $f2 ="../output/01_transfrom_kinsae_gene_to_symbol.txt";#输入的是start
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 ="../output/06_kinase_not_in_the_network.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 ="../output/06_kinase_in_the_network.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 ="../output/06_fisher_need_data.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
print $O1 "symbol\n";
print $O2 "symbol\n";
print $O3 "symbol\tkinase_in_the_top_final\tgene_not_kinase_in_the_top\tkinase_out_of_top\tgene_not_kinase_out_of_the_top\n";
my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^query/){
        my $query = $f[0];
        my $symbol = $f[1];
        $hash1{$symbol}=1;
    }
}

# my @kinases;
while(<$I2>)
{
    chomp;
    my @f= split /\t/;
    unless(/^query/){
        my $query = $f[0];
        my $symbol = $f[1];
        if (exists $hash1{$symbol}){ #在network 的kinase
        print $O2 "$symbol\n";
            $hash2{$symbol}=1;
            $hash3{$symbol}=1;
            # push @kinases,$symbol;
        }
        else{ #不在network 的kinase
            print $O1 "$symbol\n";
        }
    }
}

foreach my $symbol(sort keys %hash2){
    my @neighbours=();
    my $f3 ="../output/RWR/rwr_result_top/${symbol}.txt";#输入的是start
    open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
    while(<$I3>)
    {
        chomp;
        my @f =split/\t/;
        my $neighbour = $f[0];
        if (exists $hash3{$neighbour}){
            push @neighbours,$neighbour;
            # print "$symbol\t$neighbour\n";
        }
       
    }
     print "$symbol\t@neighbours\n";
    # my $kinase_in_the_top = @neighbours;
    # my $gene_not_kinase_in_the_top = 139 - $kinase_in_the_top; #TOP数目一共为139个。
    # my $kinase_out_of_top = 69 - $kinase_in_the_top;
    # my $gene_not_kinase_out_of_the_top = 13852 -$gene_not_kinase_in_the_top - $kinase_out_of_top;
    # my $kinase_in_the_top_final = $kinase_in_the_top -1; #因为top的第一个是自己，所以再减去一个
    # my $output = "$symbol\t$kinase_in_the_top_final\t$gene_not_kinase_in_the_top\t$kinase_out_of_top\t$gene_not_kinase_out_of_the_top";
    # print $O3 "$output\n";

}


