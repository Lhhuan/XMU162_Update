#/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz 在"cis_1MB","cis_10MB","trans_1MB","trans_10MB"
#按照相同染色体位置，与"caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"取交集得，../output/xQTL_merge/${final_type}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${final_type}_${QTL1}_${QTL2}.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
my $fo1 = "../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_1MB_eQTL_cis_trans_overlap.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_10MB_eQTL_cis_trans_overlap.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O1 "chr\tpos\tcis_1MB_emplambda\ttrans_1MB_emplambda\n";
print $O2 "chr\tpos\tcis_10MB_emplambda\ttrans_10MB_emplambda\n";

my (%hash1,%hash2,%hash3,%hash4);
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
        #------------------------------------------1MB
        if ($cis_or_trans =~/cis_1MB/){
            $hash1{$k}=$emplambda;
        }
        elsif($cis_or_trans =~/trans_1MB/){
            $hash2{$k}=$emplambda;
        }
        #-----------------------------------------10MB
        elsif($cis_or_trans =~/cis_10MB/){
            $hash3{$k}=$emplambda;
        }
        else{
            $hash4{$k}=$emplambda;
        }
    }
}

foreach my $k1(sort keys %hash1){
    if (exists $hash2{$k1}){
        my $v1 = $hash1{$k1};
        my $v2 =$hash2{$k1};
        print $O1 "$k1\t$v1\t$v2\n";
    }
}

foreach my $k2(sort keys %hash3){
    if (exists $hash4{$k2}){
        my $v1 = $hash3{$k2};
        my $v2 =$hash4{$k2};
        print $O2 "$k2\t$v1\t$v2\n";
    }
}

close($I1);
close($O1);
close($O2);