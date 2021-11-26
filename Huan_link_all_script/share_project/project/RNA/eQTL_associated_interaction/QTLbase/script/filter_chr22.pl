#将../output/01_all_kinds_QTL.txt 中chr22过滤得../output/chr22.txt,过滤chr22中eur得../output/chr22_eur.txt
#得所有的具有 pop 和tissue的文件得../output/merge_QTL_all_QTLtype_pop.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my (%hash1, %hash2, %hash3,%hash5);
my $f1 = "../data/QTLbase_download_data_sourceid.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../output/01_all_kinds_QTL.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $fo1 = "../output/chr22.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/chr22_eur.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 = "../output/merge_QTL_all_QTLtype_pop.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;

while(<$I1>)
{
    chomp;
    unless(/PMID/){
        # print "$file\n";
        my @f =split/\t/;
        my $Sourceid = $f[1];
        my $Tissue =$f[3];
        my $Population = $f[4];
        my $v = "$Tissue\t$Population";
        $hash5{$Sourceid}=$v; #用于为$f2添加$Tissue和$Population
    }
}


while(<$I2>)
{
    chomp;
    if(/SNP_chr/){
        print $O1 "$_\tTissue\tPopulation\n";
        print $O2 "$_\tTissue\tPopulation\n";
        print $O3 "$_\tTissue\tPopulation\n";
    }
    else{
        my @f =split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $Sourceid =$f[-2];
        if (exists $hash5{$Sourceid}){
            my $v = $hash5{$Sourceid};
            print $O3 "$_\t$v\n";
            if ($snp_chr =~/\b22\b/){
                print $O1 "$_\t$v\n";
                my @t =split/\t/,$v;
                my $pop = $t[1];
                if ($pop =~/EUR/){
                    print $O2 "$_\t$v\n";
                }
            }
        }
    }
}