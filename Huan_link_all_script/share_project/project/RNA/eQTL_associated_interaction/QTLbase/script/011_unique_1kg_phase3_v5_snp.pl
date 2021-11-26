#将"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz" 的 SNP 提出来得 在一起得../output/011_eur_1kg_phase3_v5_hg19_snp.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "../output/011_eur_1kg_phase3_v5_hg19_snp.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "CHROM\tPOS\tID\n";

    # print "$pop\n";
my $f1 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID.vcf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


while(<$I1>)
{
    chomp;
    unless(/^#/){
        my @f = split/\s+/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $ID =$f[2]; #snp chr
        my $output = "$CHROM\t$POS\t$ID";
        print $O1 "$output\n";
    }
}