#将../../../output/ALL_mQTL/cis_trans/NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_mQTL.txt.gz 中的QTLbase中../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz的mQTL数据过滤出来，
#得../../../output/ALL_mQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_mQTL.txt.gz
#并得汇总文件"../../../output/ALL_mQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_mQTL.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my $f1 = "../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


my (%hash1,%hash2,%hash3,%hash4);
# my %name;
while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $snp_chr =$f[0];
        my $snp_pos =$f[1];
        my $Pvalue =$f[6];
        my $QTL_type = $f[-5];
        my $cis_trans_1MB = $f[-2];
        my $cis_trans_10MB = $f[-1];
        if($QTL_type=~/\bmQTL\b/){#----filter mQTL
        # print "$QTL_type\n";
            my $k= "$snp_chr\t$snp_pos";
            if ($cis_trans_1MB =~/cis/){ #cis_1MB
                my $name1 ="cis_1MB";
                $hash1{$k}{$name1}=1;
            }
            else{   #if($cis_trans_1MB =~/trans/) #trans_1MB
                my $name ="trans_1MB";
               $hash1{$k}{$name}=1;
            }
            #----------------------------------
            if ($cis_trans_10MB =~/cis/){ ##cis_10MB
                my $name ="cis_10MB";
                $hash1{$k}{$name}=1;
            }
            else{   #if($cis_trans_10MB =~/trans/) ##trans_10MB
                my $name ="trans_10MB";
                $hash1{$k}{$name}=1;
            }
        }
    }
}


my $fo2 = "../../../output/ALL_mQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_mQTL.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "emplambda\tt\tchr\tcis_or_trans\n";
my @qtl=("cis_1MB","trans_1MB","cis_10MB","trans_10MB");
foreach my $name(@qtl){
    my $f2 = "../../../output/ALL_mQTL/cis_trans/NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_mQTL.txt.gz";
    # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    my $fo1 = "../../../output/ALL_mQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_mQTL.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;

    print "start\t$name\n";
    my $namew= $name;
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
            my $k1 = "$chr\t$pos";
            if (exists $hash1{$k1}{$namew}){
                # print $O1 "$_\t$namew\t$name\n";
                print $O1 "$_\n";
                print $O2 "$_\t$namew\n";

            }
            
        }
    }
    close($O1);
    print "finish\t$name\n";
}