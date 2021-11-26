#将"../../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz" 处理成$cis_or_trans:$emplambda的格式，并将位置处理为bed
#得../../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use 5.010;


my $f1 = "../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../../output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL_normalied.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my %hash1;
print $O1 "Chr\tstart\tend\tpos\templambda\n";
# my %name;
while(<$I1>)
{
    chomp;
    unless(/^emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $pos =$f[1];
        my $start = $pos;
        my $end = $start +1; #0-based single pos to bed
        my $chr =$f[2];
        my $final_chr = "chr${chr}";
        my $cis_or_trans = $f[3];
        my $k = "$final_chr\t$start\t$end\t$pos\t$emplambda";
        print $O1 "$k\n";
    }
}