#将../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 按照特定距离将emplambda取平均
#得../output/average_emplambda_xQTL.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.chr22.vcf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/share/data0/QTLbase/data/eQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
#------------------------------------

my $fo1 = "../output/QTLbase_chr22_region.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "../output/QTLbase_chr22_and_1kg.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

my (%hash1,%hash2);


while(<$I2>)
{
    chomp;
    if(/SNP_chr/){
        print $O1 "$_\n";
        print $O2 "$_\n";
    }
    else{
        my @f = split/\t/;
        my $SNP_chr =$f[0];
        my $POS =$f[1]; 
        if ($POS >16000000 && $POS<17000000){ #20000个
            print $O1 "$_\n";
            print $O2 "$_\n";
            $hash1{$POS}=1;


        }
    
    }
}

while(<$I1>)
{
    chomp;
    unless(/^#/){
        my @f = split/\s+/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $ID =$f[2]; #snp chr
        print "$POS\n";
        if ($POS >16000000 && $POS<17000000){
            unless(exists $hash1{$POS}){
                print $O2 "$CHROM\t$POS\t\t\t\t\t0.05\t1kg\n";
            }
        }
    }
}