 # 用../output/all_1kg_phase3_v5_hg19_snp.txt.gz 补全../output/merge_QTL_all_QTLtype_pop.txt.gz 得../output/merge_QTL_all_QTLtype_pop_1kg_Completion.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "../output/merge_QTL_all_QTLtype_pop_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

my $f1 = "../output/merge_QTL_all_QTLtype_pop.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $f2 = "../output/all_1kg_phase3_v5_hg19_snp.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

print $O1 "SNP_chr\tSNP_pos\tQTL_type\tPvalue\n";
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $SNP_chr =$f[0];
        my $SNP_pos =$f[1];
        my $Pvalue =$f[6];
        my $QTL_type = $f[-3];
        my $k = "$SNP_chr\t$SNP_pos\t$QTL_type";
        $hash1{$k}=1;
        print $O1 "$k\t$Pvalue\n";
    }
}


my @qtl_type = ("caQTL","cerQTL","eQTL","edQTL","hQTL","lncRNAQTL","mQTL","metaQTL","miQTL","pQTL","reQTL","riboQTL","sQTL");
while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        foreach my $QTL (@qtl_type) {
            my $k = "$CHROM\t$POS\t$QTL";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
            }
        }

    }
}



