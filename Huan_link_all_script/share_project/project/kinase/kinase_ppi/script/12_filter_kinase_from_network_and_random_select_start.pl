#将../output/03_final_all_network_gene_symbol.txt 中的../output/08_kinase_in_the_network.txt
#过滤掉，并随机选择69个基因。得../output/12_random_start.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use File::Copy;

my $f1 ="../output/08_kinase_in_the_network.txt";
my $f2 ="../output/03_final_all_network_gene_symbol.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 ="../output/12_random_start.txt"; 
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^symbol/){
        my $symbol = $f[0];
        $hash1{$symbol}=1;
    }
}

my @genes=();
while(<$I2>)
{
    chomp;
    my @f= split /\t/;
    unless(/^Gene/){
        # my $query = $f[0];
        my $symbol = $f[1];
        unless(exists $hash1{$symbol}){
            push @genes,$symbol;
        }
    }
}
# my @array = @genes;
my $num = 69;
my @new = randomElem ( $num, @genes ) ; # pick any $num from @array ，把$num和@array传递给子程序。这里是用的值传递。还有一种方式是引用传递，相当于硬链接
my $out = join ("\n",@new);
print $O1 "$out\n";

sub randomElem { #随机取
    my ($want, @array) = @_ ;
    my (%seen, @ret);
    while ( @ret != $want ) {
    my $num = abs(int(rand(@array))); #@array 是指数组的长度，而$#array是指最后一个索引，由于rand的特殊性，如果用$#array会导致取不到最后一个值。
        if ( ! $seen{$num} ) { 
        ++$seen{$num};
        push @ret, $array[$num];
        }
    }
    return @ret;     
}
