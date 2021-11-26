#将../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz 中的QTLbase中../output/merge_QTL_all_QTLtype_pop.txt.gz的数据过滤出来，
#得../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;


my (%hash1,%hash2,%hash4);
my $f1 = "/share/data0/QTLbase/data/eQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../../output/ALL_eQTL/QTLbase_all_eQTL.bed";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

while(<$I1>)
{
    chomp;
    unless(/SNP_chr/){
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $Pvalue =$f[-2];
        if ($Pvalue <5e-8){
            # if ($snp_chr  eq 1 $snp_chr  eq 1){
                my $start = $snp_pos;
                my $end = $start +1;
                my $value =1;
                my $output ="$snp_chr\t$start\t$end\t$value";
                unless(exists $hash1{$output}){
                    $hash1{$output} =1;
                    print $O1 "$output\n";
                }
            # }
        }
    }
}

        