#将../output/ALL_${xQTL}/NHPoisson_emplambda_interval_${interval}cutoff_${cutoff}_all_${xQTL}.txt.gz 中的QTLbase中../output/merge_QTL_all_QTLtype_pop.txt.gz的数据过滤出来，
#得../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${interval}_cutoff_${cutoff}_all_${xQTL}.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my (%hash1,%hash2,%hash4);
my $f1 = "../output/merge_QTL_all_QTLtype_pop.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


while(<$I1>)
{
    chomp;
    unless(/SNP_chr/){
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $QTL_type = $f[8];
        my $k1 = "$snp_chr\t$snp_pos\t$QTL_type";
        my $k2 = "$snp_chr\t$snp_pos";
        $hash1{$k1}=1;
        $hash2{$k2}=1;

        
    }
}


my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL","cerQTL","lncRNAQTL","metaQTL","miQTL","riboQTL");
my @cutoffs=(7.3,7);
my @intervals =(1000,2000,3000,4000,5000);
# my $xQTL = "eQTL" ;
foreach my $xQTL(@QTLs){
    foreach my $cutoff(@cutoffs){
        foreach my $interval(@intervals){
            my $f2 = "../output/ALL_${xQTL}/NHPoisson_emplambda_interval_${interval}cutoff_${cutoff}_all_${xQTL}.txt.gz";
            # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            my $fo1 = "../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${interval}_cutoff_${cutoff}_all_${xQTL}.txt.gz";
            open my $O1, "| gzip >$fo1" or die $!;
            print "start\t$xQTL\n";
            while(<$I2>)
            {
                chomp;
                if(/emplambda/){
                    print $O1 "$_\n";
                }
                else{
                    my @f = split/\t/;
                    my $emplambda =$f[0];
                    my $t =$f[1]; 
                    my $pos = $t; #snp pos
                    my $chr =$f[2]; #snp chr
                    my $k1 = "$chr\t$pos\t$xQTL";
                    my $k2 = "$chr\t$pos";
                    if ($xQTL =~/\bQTL\b/){ #---all QTL
                        if(exists $hash2{$k2}){
                            print $O1 "$_\n";
                        }
                    }
                    else{ #------other QTL
                        if (exists $hash1{$k1}){
                            print $O1 "$_\n";
                        }
                    }
                    
                }
            }
            print "finish\t$xQTL\n";
        }
    }
}