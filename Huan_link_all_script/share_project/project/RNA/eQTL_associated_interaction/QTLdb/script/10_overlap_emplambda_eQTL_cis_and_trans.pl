#/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz 在"cis_1MB","cis_10MB","trans_1MB","trans_10MB"
#按照相同染色体位置，与"caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"取交集得，../output/xQTL_merge/${final_type}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${final_type}_${QTL1}_${QTL2}.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $t =$f[1]; 
        my $pos = $t; #snp pos
        my $chr =$f[2]; #snp chr
        my $cis_or_trans =$f[3];
        my $k = "$chr\t$pos";
        unless ($emplambda =~/NA/){
            $hash1{$k}{$cis_or_trans}=$emplambda;
        }

    }
}




# # my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
# my @QTL1s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");
my @QTL2s = ("caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");
my @types = ("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my $QTL1 = "eQTL";
foreach my $QTL2(@QTL2s){
    # print "$QTL1\t$QTL2\n";
    foreach my $type(@types){
        print "$QTL1\t$QTL2\t$type\n";
        my $final_type = ucfirst($type);
        my $f2 = "../output/ALL_${QTL2}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL2}.txt.gz";
        # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
        open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
        my $fo1 = "../output/xQTL_merge/${final_type}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${final_type}_${QTL1}_${QTL2}.txt.gz";
        open my $O1, "| gzip >$fo1" or die $!;
        print $O1 "${QTL1}_emplambda\t${QTL2}_emplambda\tpos\tchr\n";
        # print "$cp_format1\n";
        
        while(<$I2>)
        {
            chomp;
            unless(/emplambda/){
                my @f = split/\t/;
                my $emplambda =$f[0];
                my $t =$f[1]; 
                my $pos = $t; #snp pos
                my $chr =$f[2]; #snp chr
                my $k = "$chr\t$pos";
                if (exists $hash1{$k}{$type}){
                    my $emplambda1 = $hash1{$k}{$type};
                    print $O1 "$emplambda1\t$emplambda\t$chr\t$pos\n";
                }
            }
        }
        close($I1);
        close($O1);
    }
}