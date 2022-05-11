#提取gene 名字得 "01_unique_gene.csv.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "goa_human.gaf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "01_unique_gene.csv.gz"; #
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "gene\n";
# print $O2 "query_gene\tensembl\n";




my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^\!/){
        my $gene=$f[2];
        unless(exists $hash1{$gene}){
            $hash1{$gene}=1;
            print $O1 "$gene\n";
        }
    }
}
