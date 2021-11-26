 #("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt按照相同染色体位置，两两取交集，
#得../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_${QTL1}_${QTL2}.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my %hash1;


# my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
my @QTL1s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");
my @QTL2s = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");

foreach my $QTL1(@QTL1s){
    foreach my $QTL2(@QTL2s){
        unless ($QTL1 =~/$QTL2/){#防止自己与自己取交集
            my $cp_format1 = "$QTL1\t$QTL2";
            my $cp_format2 = "$QTL2\t$QTL1";
            unless(exists $hash1{$cp_format1} || exists $hash1{$cp_format2} ){
                my $f1 = "../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL1}.txt";
                open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
                my $f2 = "../output/ALL_${QTL2}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${QTL2}.txt";
                open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
                my $fo1 = "../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz";
                open my $O1, "| gzip >$fo1" or die $!;
                print $O1 "${QTL1}_emplambda\t${QTL2}_emplambda\tpos\tchr\n";
                print "$cp_format1\n";
                my %hash2;
                while(<$I1>)
                {
                    chomp;
                    unless(/emplambda/){
                        my @f = split/\t/;
                        my $emplambda =$f[0];
                        my $t =$f[1]; 
                        my $pos = $t; #snp pos
                        my $chr =$f[2]; #snp chr
                        my $k = "$chr\t$pos";
                        unless ($emplambda =~/NA/){ 
                            $hash2{$k}=$emplambda;

                        }
                    }
                }
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
                        unless ($emplambda =~/NA/){ 
                            if(exists $hash2{$k} ){
                                my$emplambda1  = $hash2{$k};
                                print $O1 "$emplambda1\t$_\n";
                            }
                        }
                    }
                }
                close($O1);
                close($I1);
                close($I2);
            $hash1{$cp_format1}=1; #cp在上面没出现过
            $hash1{$cp_format2}=1;
            }
        }
    }
}