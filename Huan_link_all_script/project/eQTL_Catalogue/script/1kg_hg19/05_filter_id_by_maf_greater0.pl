#将"./ID_MAF/EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.maf.id.txt.gz"  中的EUR提出来，得EUR_sample_list.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "EUR_id_maf_greater0.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my (%hash1,%hash2);

print $O1 "chr_pos_ref_alt\tmaf\n";


my $f2 = "/home/huanhuan/project/eQTL_Catalogue/script/1kg_hg19/1kg.phase3.v5.shapeit2.eur.hg19.all.maf.id.vcf.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f =split/\t/;
        my $id =$f[0];
        my $maf =$f[1];
        # print "$_\n";
        if($maf =~/\,/){
            my @ts =split/,/,$maf;
            foreach my $t(@ts){
                unless(exists $hash1{$_}) {
                    $hash1{$_}=1;
                    print $O1 "$_\n";
                }        
            }
        }
        else{
            if($maf >0){ 
                unless(exists $hash1{$_}) {
                    $hash1{$_}=1;
                    print $O1 "$_\n";
                }
            }
        }
    }
}

