#将/home/chaoqun/phase3/integrated_call_samples_v3.20130502.ALL.panel 中的EUR提出来，得EUR_sample_list.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "/share/Projects/huanhuan/ref_data/dbSNP/hg38/All_20180418.vcf.gz";
# my $f1 = "12345.vcf.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "./ID_MAF/snp/EUR_snp_maf_more_than0.05.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my (%hash1,%hash2);

print $O1 "chr:pos:ref:alt\tmaf\trsid\n";
while(<$I1>)
{
    chomp;
    unless(/^#/){
        my @f =split/\t/;
        my $CHROM =$f[0];
        my $POS = $f[1];
        my $rs_ID =$f[2];
        my $REF =$f[3];
        my $ALT =$f[4];
        my $id = "$CHROM:$POS:$REF:$ALT";
        $hash1{$id}=$rs_ID;
        # print "$id\t$rs_ID\n";
    }
}

my @numbers=(1..22);
foreach my $i(@numbers){
    print "$i\n";
    my $f2 = "./ID_MAF/EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.maf.id.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f =split/\t/;
            my $id =$f[0];
            my $maf =$f[1];
            # print "$_\n";
            if($maf >0.05){                
                if(exists $hash1{$id}){
                    my $rsid =$hash1{$id};
                    print $O1 "$_\t$rsid\n";
                }
            }
        }
    }
}

