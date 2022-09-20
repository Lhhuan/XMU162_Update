
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

my $f1 = "../output/huan_GATNE/train.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../output/huan_GATNE/test.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $f3 = "../output/huan_GATNE/valid.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";

# my $fo1 = "../output/02_hotspot_egene_idx.txt.gz"; #
# open my $O1, "| gzip >$fo1" or die $!;

my (%hash1,%hash2,%hash3);
while(<$I1>)
{
    chomp;
    my @f=split/\s+/;
    my $k1=join("\t",@f[0..2]);
    $hash1{$k1}=1;
    # print "$k1\n";
}

while(<$I2>)
{
    chomp;
    my @f=split/\s+/;
    my $k2=join("\t",@f[0..2]);
    $hash2{$k2}=1;
    if(exists $hash1{$k2}){
        print "$k2\ttest\n";        
    }
    # print "$k2\n";
}

while(<$I3>)
{
    chomp;
    my @f=split/\s+/;
    my $k3=join("\t",@f[0..2]);
    $hash3{$k3}=1;
    if(exists $hash2{$k3}){
        print "$k3\n";
    }
}