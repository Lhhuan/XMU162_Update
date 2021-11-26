#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/5mer/5communities.csv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/5mer/5communities.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my %hash1;

while(<$I1>)
{
    chomp;
    if(/^>/){
        my @f = split/,/;
        my $p = $f[0];
        my $class =$f[1];
        $p=~s/>//g;
        # print "$p\n";
        my @ps= split/:/,$p;
        my $chr = $ps[0];
        my @ses = split/-/,$ps[1];
        my $start = $ses[0];
        my $end = $ses[1];
        print $O1 "$chr\t$start\t$end\t$class\n";

    }
}

close($I1);
close($O1);

my @factors = ("TFBS","CHROMATIN_Accessibility","CTCF","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3");
my $file_name = "Whole_Blood_segment_hotspot_cutoff_0.176.bed.gz";
my $factor_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/0.176";


foreach my $factor(@factors){
    my $factor_file = "$factor_dir/${factor}_${file_name}";
    my $command = "bedtools intersect -a $fo1  -b $factor_file -wa |sort -u |gzip >/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/cluster5/5mer/factor_anno/${factor}.bed.gz" ;
    print "$factor_file\n";
    system $command;
}
