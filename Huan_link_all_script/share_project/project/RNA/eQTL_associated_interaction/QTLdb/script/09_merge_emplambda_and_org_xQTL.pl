#将../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz 和 ../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz merge在一起得
#cis和trans 信息文件为../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz
#得snp和gene 位于同一染色体上的 cis, trans 及distance 信息为 ../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my (%hash1,%hash3,%hash4);
my $f1 = "../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $f2 = "../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

while(<$I1>)
{
    chomp;
    unless(/SNP_chr/){
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $Mapped_gene =$f[2];
        my $trait_chr = $f[3];
        my $start = $f[4];
        my $end  = $f[5];
        my $Pvalue =$f[6];
        my $QTL_type = $f[8];
        my $cis_trans_1MB = $f[-2];
        my $cis_trans_10MB =$f[-1];
        my $k1 = "$snp_chr\t$snp_pos\t$QTL_type";
        my $v1 = $Pvalue ;
        push @{$hash1{$k1}},$v1;
    }
}

#-------------------------------------select min pvalue
foreach my $k1 (sort keys %hash1){
    my @vs = @{$hash1{$k1}};
    my $length = @vs;
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs;
    if ($length >1){
        my $min_p = min @vs;
        my $k3 = "$k1\t$min_p";
        $hash3{$k3}=1;
    }
    else{
        my $k3 = "$k1\t$vs[0]";
        $hash3{$k3}=1;
    }
}
#----------------------------------------------

while(<$I2>)
{
    chomp;
    unless(/SNP_chr/){
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $Mapped_gene =$f[2];
        my $trait_chr = $f[3];
        my $start = $f[4];
        my $end  = $f[5];
        my $Pvalue =$f[6];
        my $QTL_type = $f[8];
        my $cis_trans_1MB = $f[-2];
        my $cis_trans_10MB =$f[-1];
        #-------------------------------------------# select the value of the min p
        my $k3 = "$snp_chr\t$snp_pos\t$QTL_type\t$Pvalue";
        if (exists $hash3{$k3}){
            my $k = "$snp_chr\t$snp_pos\t$QTL_type";
            my $v = "$Mapped_gene\t$trait_chr\t$start\t$end\t$Pvalue\t$cis_trans_1MB\t$cis_trans_10MB";
            push @{$hash4{$k}},$v;
        }
    }
}

# my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL");

# my $xQTL = "eQTL" ;
foreach my $xQTL(@QTLs){
    my $f3 = "../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz";
    # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
    my $fo1 = "../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "emplambda\tt\tchr\tMapped_gene\ttrait_chr\tstart\tend\tPvalue\tcis_trans_1MB\tcis_trans_10MB\n";
    my $fo2 = "../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}_cis_trans_distance.txt.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    print $O2 "emplambda\tt\tchr\tMapped_gene\ttrait_chr\tstart\tend\tPvalue\tcis_trans_1MB\tcis_trans_10MB\tdistance\n";
    while(<$I3>)
    {
        chomp;
        unless(/emplambda/){
            my @f = split/\t/;
            my $emplambda =$f[0];
            my $t =$f[1]; 
            my $pos = $t; #snp pos
            my $chr =$f[2]; #snp chr
            my $k = "$chr\t$pos\t$xQTL";
            unless ($emplambda =~/NA/){        
                if (exists $hash4{$k}){
                    my @vs = @{$hash4{$k}};
                    my %hash5;
                    @vs = grep { ++$hash5{$_} < 2 } @vs;
                    foreach my $v (@vs){
                        my $output1 = "$_\t$v";
                        print $O1 "$output1\n";
                        my @t =split/\t/,$v;
                        my $Mapped_gene =$t[0];
                        my $trait_chr =$t[1];
                        my $start =$t[2];
                        my $end =$t[3];
                        my $Pvalue = $t[4];
                        my $cis_trans_1MB =$t[5];
                        my $cis_trans_10MB =$t[6];
                        if ($chr =~/$trait_chr/){ #snp in gene
                            if ($start>=$pos && $end<= $pos){ #SNP in gene 
                                print $O2 "$output1\t0\n";
                            }
                            else{
                                my $length_start = abs($start -$pos);
                                my $length_end = abs($end-$pos);
                                if ($length_start <=$length_end){  #select the min distance of gene and snp
                                    print $O2 "$output1\t$length_start\n";
                                }
                                else{
                                    print $O2 "$output1\t$length_end\n";
                                }
                            }
                            
                        }
                    }
                }
            }
        }
    }
}