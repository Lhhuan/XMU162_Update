#"/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz"  注释 "03_entrezgene_symbol_ensembl.txt.gz" 的位置，得"03_entrezgene_symbol_ensembl_pos.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1= "/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2= "03_entrezgene_symbol_ensembl.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); 
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";

my $fo1 = "03_entrezgene_symbol_ensembl_pos.txt.gz"; #
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "entrezgene\tensembl\tsymbol\tChr\tstart\tend\n";

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/Gene/){
        my @f=split/\t/;
        my $ensg =$f[0];
        my $start =$f[1];
        my $end=$f[2];
        my $chr=$f[3];
        unless($chr=~/CHR|GL|KI|MT|X|Y/){
            my $v= "chr${chr}\t$start\t$end";
            $hash1{$ensg}=$v;
        }
    }
}

while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    unless(/^entrezgene/){
        my $gene =$f[0];
        my $ensg = $f[1];
        if(exists $hash1{$ensg}){
            my $pos = $hash1{$ensg};
            print $O1 "$_\t$pos\n";
        }
    }
}
