#("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL"#中，用../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz按照相同染色体位置，两两取交集，
#得 QTL1为trans, QTL2为cis文件 ../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_trans_${QTL1}_cis_${QTL2}.txt.gz
#和../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_trans_${QTL1}_cis_${QTL2}.txt.gz

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
                my $f1 = "../output/ALL_${QTL1}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL1}_cis_trans.txt.gz";
                # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
                open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
                my $f2 = "../output/ALL_${QTL2}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${QTL2}_cis_trans.txt.gz";
                # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
                open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                my $fo1 = "../output/xQTL_merge/Cis_Trans_1MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_trans_${QTL1}_cis_${QTL2}.txt.gz";
                open my $O1, "| gzip >$fo1" or die $!;
                print $O1 "trans_${QTL1}_emplambda\tcis_${QTL2}_emplambda\tpos\tchr\n";
                my $fo2 = "../output/xQTL_merge/Cis_Trans_10MB/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_trans_${QTL1}_cis_${QTL2}.txt.gz";
                open my $O2, "| gzip >$fo2" or die $!;
                print $O2 "trans_${QTL1}_emplambda\tcis_${QTL2}_emplambda\tpos\tchr\n";
                print "$cp_format1\n";
                my %hash2;
                my %hash3;
                while(<$I1>)
                {
                    chomp;
                    unless(/emplambda/){
                        my @f = split/\t/;
                        my $emplambda =$f[0];
                        my $t =$f[1]; 
                        my $pos = $t; #snp pos
                        my $chr =$f[2]; #snp chr
                        my $cis_trans_1MB =$f[-2];
                        my $cis_trans_10MB =$f[-1];
                        my $k = "$chr\t$pos";
                        unless ($emplambda =~/NA/){
                #-----------------------------------------------
                            if ($cis_trans_1MB =~/trans/){
                                    $hash2{$k}=$emplambda;
                            }
                            
                    #-----------------------------------------------
                            if ($cis_trans_10MB =~/trans/){ 
                                    $hash3{$k}=$emplambda;
                            }
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
                        my $cis_trans_1MB =$f[-2];
                        my $cis_trans_10MB =$f[-1];
                        my $k = "$chr\t$pos";
                        unless ($emplambda =~/NA/){ 
                            #----------------------------------1MB as cis
                            if ($cis_trans_1MB =~/cis/){
                                if(exists $hash2{$k} ){
                                    my$emplambda1  = $hash2{$k};
                                    print $O1 "$emplambda1\t$emplambda\t$pos\t$chr\n";
                                }
                            }
                            #-----------------------------------------10MB as cis
                            if($cis_trans_10MB =~ /cis/){
                                if(exists $hash3{$k} ){
                                    my$emplambda1  = $hash3{$k};
                                    print $O2 "$emplambda1\t$emplambda\t$pos\t$chr\n";
                                }                                
                            }
                        }
                    }
                }
                close($O1);
                close($O2);
                close($I1);
                close($I2);
            $hash1{$cp_format1}=1; #cp在上面没出现过
            $hash1{$cp_format2}=1;
            }
        }
    }
}