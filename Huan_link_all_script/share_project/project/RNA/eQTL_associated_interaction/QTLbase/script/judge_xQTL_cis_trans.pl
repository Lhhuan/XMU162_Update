#将../output/merge_QTL_all_QTLtype_pop.txt.gz 按照1MB和10MB划分cis和trans， 分别得../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my (%hash1, %hash2, %hash3,%hash5);
my $f1 = "../output/merge_QTL_all_QTLtype_pop.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;


while(<$I1>)
{
    chomp;
    if(/SNP_chr/){
        print $O1 "$_\tcis_trans_1MB\tcis_trans_10MB\n";
    }
    else{
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $trait_chr = $f[3];
        my $start = $f[4];
        my $end  = $f[5];
        my $spos = $snp_pos;
        my $up_snp_1MB =$spos-1000000;
        my $down_snp_1MB = $spos+1000000;
        #--------------
        my $up_snp_10MB =$spos-10000000;
        my $down_snp_10MB = $spos+10000000;
        #------------------------------1MB 
        print $O1 "$_\t";
        if ($snp_chr == $trait_chr  && $start>=$up_snp_1MB && $end<= $down_snp_1MB ){
            print $O1 "cis\t";
        }
        else{
            print $O1 "trans\t";
        }
        #--------------------------10MB
        if ($snp_chr == $trait_chr  && $start>=$up_snp_10MB && $end<= $down_snp_10MB ){
            print $O1 "cis\n";
        }
        else{
            print $O1 "trans\n";
        }       
    }
}

