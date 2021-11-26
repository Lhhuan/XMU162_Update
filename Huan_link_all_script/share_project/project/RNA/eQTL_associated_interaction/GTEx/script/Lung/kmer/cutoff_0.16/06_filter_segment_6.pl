#过滤../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz中6-5000 bp的片段,得../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6_5000/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18/Lung_segment_hotspot_cutoff_0.16.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#----------------------
my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/Lung_segment_hotspot_cutoff_0.16.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;



while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $SNP_chr =$f[0];
    my $start =$f[1];
    my $end =$f[2];
    my $length = $end-$start;
    if($length >5){
        print $O1 "$_\n";
    }
}

