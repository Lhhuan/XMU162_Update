#利用../output/0011_QTLbase_download_data_sourceid_QC.txt，提取/share/data0/QTLbase/data/xQTL, 得/share/data0/QTLbase/2020_12_2_QC/xQTLd_QC.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;




my @xQTLs = ("caQTL.txt.gz","cerQTL.txt.gz","edQTL.txt.gz","eQTL.txt.gz","hQTL.txt.gz","lncRNAQTL.txt.gz","metaQTL.txt.gz",
    "miQTL.txt.gz","mQTL.txt.gz","pQTL.txt.gz","reQTL.txt.gz","riboQTL.txt.gz","sQTL.txt.gz");

my $f1 = "../output/0011_QTLbase_download_data_sourceid_QC.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^PMID/){
        my @f = split/\t/;
        my $PMID =$f[0];
        my $Sourceid = $f[1];
        $hash1{$Sourceid}=1;
    }
}

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $pm = Parallel::ForkManager->new(15);
foreach my $xQTL(@xQTLs){
    my $pid = $pm->start and next; #开始多线程
    my $fo1 = "/share/data0/QTLbase/2020_12_2_QC/${xQTL}";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;

    my $f2 = "/share/data0/QTLbase/data/${xQTL}";
    # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

    while(<$I2>)
    {
        chomp;
        
        if(/^SNP_chr/){
            print $O1 "$_\n";
        }
        else{
            my @f = split/\t/;
            my $Sourceid =$f[-1];
            if(exists $hash1{$Sourceid}){
                print $O1 "$_\n";
            }
        }
    }
    $pm->finish;  #多线程结束
}