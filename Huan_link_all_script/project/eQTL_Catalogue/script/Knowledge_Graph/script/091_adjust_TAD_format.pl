#Adjust "../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz"得 "../output/edges_annotation/09_Adjust_success_egene_hotspot_TAD.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my (%hash1,%hash2);

my $fo1 = "../output/edges_annotation/09_Adjust_success_egene_hotspot_TAD.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Hotspot_gene\tTAD\n";
my $f1= "../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

while(<$I1>)
{
    chomp;
    my @f=split/\t/;
    my $gene_p =join("\t",@f[0..2]);
    my $hotspot_p =join("_",@f[3..5]);
    my $egene =$f[6];
    my $TAD = join("_",@f[7..9]);
    my $tissue =$f[10];
    my $k="$hotspot_p:$egene";
    my $v= "${TAD}:${tissue}";
    push @{$hash1{$k}},$v;
}


foreach my $k(sort keys %hash1){
    my @vs=@{$hash1{$k}};
    @vs=uniq(@vs);
    my $v=join(";",@vs);
    print $O1 "$k\t$v\n";

}