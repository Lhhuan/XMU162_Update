#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;
my $f1 = "/share/data0/QTLbase/huan/hic/4DNFI2TK7L2F_5000.ginteractions.tsv.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 

my $fo2 = "/share/data0/QTLbase/huan/hic/4DNFI2TK7L2F_5000.ginteractions_cross_chr_adjust.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $chr1 = $f[0];
    my $start1 = $f[1];
    my $end1 = $f[2];
    my $chr2= $f[3];
    my $start2= $f[4];
    my $end2=$f[5];
    my $signal=$f[6];
    # unless($chr1 =~/X|Y/){
        # if($chr1 eq $chr2){
            $chr1="chr${chr1}";
            $chr2="chr${chr2}";
            my $output1= "$chr1\t$start1\t$end1\t$signal";
            my $output2= "$chr2\t$start2\t$end2\t$signal";
            unless(exists $hash1{$output1}){
                $hash1{$output1}=1;
                print $O2 "$output1\n";
            }
            unless(exists $hash1{$output2}){
                $hash1{$output2}=1;
                print $O2 "$output2\n";
            }
        # }
    # }

}

close($O2);
system "zless $fo2 |sort -k1,1 -k2,2n |gzip > /share/data0/QTLbase/huan/hic/4DNFI2TK7L2F_5000.ginteractions_cross_chr_adjust_sorted.bed.gz ";







