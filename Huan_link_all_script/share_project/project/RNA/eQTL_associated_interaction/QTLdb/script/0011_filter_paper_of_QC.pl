# 筛选../output/unique_journal_source_IF.txt中IF>=5的journal,从"/share/data0/QTLbase/data/QTLbase_download_data_sourceid.txt"筛选得../output/0011_QTLbase_download_data_sourceid_QC.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "../output/0011_QTLbase_download_data_sourceid_QC.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
my $f1 = "../output/unique_journal_source_IF.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $f2 = "/share/data0/QTLbase/data/QTLbase_download_data_sourceid.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^Journal/){
        my @f = split/\t/;
        my $Journal =$f[0];
        my $IF = $f[1];
        if($IF >=5){
            $hash1{$Journal}=1;
        }
        
    }
}


while(<$I2>)
{
    chomp;
    
    if(/^PMID/){
        print $O1 "$_\n";
    }
    else{
        my @f = split/\t/;
        my $Journal =$f[6];
        if(exists $hash1{$Journal}){
            print $O1 "$_\n";
        }
    }
}