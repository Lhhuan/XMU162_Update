 # 用../../../output/all_1kg_phase3_v5_hg19_snp.txt.gz 补全../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz 得/home/huanhuan/project/RNA/cerQTL_associated_interaction/QTLbase/output/ALL_cerQTL/cis_trans/cis_1MB_cerQTL_pop_1kg_Completion.txt.gz
#trans_1MB_cerQTL_pop_1kg_Completion.txt.gz, trans_1MB_cerQTL_pop_1kg_Completion.txt.gz, cis_10MB_cerQTL_pop_1kg_Completion.txt.gz,
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $f1 = "../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $f2 = "../../../output/all_1kg_phase3_v5_hg19_snp.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_cerQTL/cis_trans/cis_1MB_cerQTL_pop_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_cerQTL/cis_trans/trans_1MB_cerQTL_pop_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O2, "| gzip >$fo2" or die $!;
my $fo3 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_cerQTL/cis_trans/cis_10MB_cerQTL_pop_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O3, "| gzip >$fo3" or die $!;
my $fo4 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_cerQTL/cis_trans/trans_10MB_cerQTL_pop_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O4, "| gzip >$fo4" or die $!;


print $O1 "SNP_chr\tSNP_pos\tQTL_type\tPvalue\n";
print $O2 "SNP_chr\tSNP_pos\tQTL_type\tPvalue\n";
print $O3 "SNP_chr\tSNP_pos\tQTL_type\tPvalue\n";
print $O4 "SNP_chr\tSNP_pos\tQTL_type\tPvalue\n";
my (%hash1,%hash2,%hash3,%hash4);
while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $SNP_chr =$f[0];
        my $SNP_pos =$f[1];
        my $Pvalue =$f[6];
        my $QTL_type = $f[-5];
        my $cis_trans_1MB = $f[-2];
        my $cis_trans_10MB = $f[-1];
        if($QTL_type=~/\bcerQTL\b/){#----filter cerQTL
        print "$QTL_type\n";
            my $k = "$SNP_chr\t$SNP_pos\t$QTL_type";
            if ($cis_trans_1MB =~/cis/){ #cis_1MB
                $hash1{$k}=1;
                print $O1 "$k\t$Pvalue\n";
            }
            else{   #if($cis_trans_1MB =~/trans/) #trans_1MB
                $hash2{$k}=1;
                print $O2 "$k\t$Pvalue\n";
            }
            #----------------------------------
            if ($cis_trans_10MB =~/cis/){ ##cis_10MB
                $hash3{$k}=1;
                print $O3 "$k\t$Pvalue\n";
            }
            else{   #if($cis_trans_10MB =~/trans/) ##trans_10MB
                $hash4{$k}=1;
                print $O4 "$k\t$Pvalue\n";
            }
        }
    }
}


while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS\tcerQTL";
#----------------------------
        unless (exists $hash1{$k}){
            print $O1 "$k\t$pvalue\n";
        }
        #-----------------
        unless (exists $hash2{$k}){
            print $O2 "$k\t$pvalue\n";
        }
        #-------------------
        unless (exists $hash3{$k}){
            print $O3 "$k\t$pvalue\n";
        }
        #-----------------
        unless (exists $hash4{$k}){
            print $O4 "$k\t$pvalue\n";
        }
    }
}



