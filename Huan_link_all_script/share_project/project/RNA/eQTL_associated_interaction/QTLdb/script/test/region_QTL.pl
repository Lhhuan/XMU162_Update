#统计../output/all_qtl_clump_locus.txt.gz 中每特定距离中QTL(tag snp的数量)，得../output/count_per_500bp_QTL_number.txt, ../output/count_per_750bp_QTL_number.txt, ../output/count_per_1000bp_QTL_number.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my (%hash1, %hash2, %hash3);

my $f1 = "../output/all_qtl_clump_locus_r_square0.5.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/EUR_Blood_mQTL.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $header = "chr\tposition\tTag_snp";
print $O1 "$header\n";


while(<$I1>)
{
    chomp;
    unless(/PMID/){
        # print "$file\n";
        my @f =split/\t/;
        my $Population =$f[-1];
        my $Tissue =$f[-2];
        my $QTL_type =$f[-3];
        my $Locus_region =$f[-5];
        my $Tag_snp =$f[-6];
        if ($Population =~/EUR/ && $Tissue =~ /\bBlood\b/ && $QTL_type =~/mQTL/){
            my @t = split/\_/,$Tag_snp;
            my $chr = $t[0];
            my $pos = $t[1];
            my $output = "$chr\t$pos\t$Tag_snp";
            unless (exists $hash1{$output}){
                $hash1{$output} =1;
                print $O1 "$output\n";
            }
        }
    }
}