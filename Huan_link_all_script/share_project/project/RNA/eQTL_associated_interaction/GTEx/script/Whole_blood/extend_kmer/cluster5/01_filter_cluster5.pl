#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group ="hotspot";
my $tissue = "Whole_Blood";


my $f1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/communities.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# my $fo1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster0/cluster0.bed.gz"; 
my $fo2 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/cluster5/cluster5.bed.gz"; 
 


# open my $O1, "| gzip >$fo1" or die $!;
open my $O2, "| gzip >$fo2" or die $!;

# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $SNP_chr =$f[0];
    my $start =$f[1];
    my $end =$f[2];
    my $cluster = $f[3];
    if($cluster ==5){
        print $O2 "$_\n";
    }
    # elsif($cluster ==5){
    #     print $O2 "$_\n";
    # }
}
# print "$group\t$cutoff\n";
      