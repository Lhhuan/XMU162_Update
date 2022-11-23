#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;


# my @sample_name=("4DNFIQYQWPF5","4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K");
# foreach my $sample(@sample_name){
my $sample = "4DNFI18Q799K";
    my $f1 = "/share/data0/QTLbase/huan/hic/${sample}_5000.ginteractions.tsv.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 

    my $fo2 = "/share/data0/QTLbase/huan/hic/${sample}_5000.ginteractions_inter_chr_adjust.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    my %hash1;

    while(<$I1>)
    {
        chomp;
        my @f= split/\t/;
        my $line_num=$.;
        my $chr1 = $f[0];
        my $start1 = $f[1];
        my $end1 = $f[2];
        my $chr2= $f[3];
        my $start2= $f[4];
        my $end2=$f[5];
        my $signal=$f[6];
        unless($chr1 =~/X|Y/){
            unless($chr2 =~/X|Y/){
                $chr1="chr${chr1}";
                $chr2="chr${chr2}";
                my $output1= "$chr1\t$start1\t$end1\t$signal\t$line_num";
                my $output2= "$chr2\t$start2\t$end2\t$signal\t$line_num";
                if($output1 eq $output2){
                    print $O2 "$output1\n";
                }
                else{
                    print $O2 "$output1\n";
                    print $O2 "$output2\n";
                }
            }
        }

    }
    close($O2);
    # my $sort_fo2= "/share/data0/QTLbase/huan/hic/${sample}_5000.ginteractions_inter_chr_adjust.bed_sorted.gz";
    # system "zless $fo2 |sort -k1,1 -k2,2n |gzip >$sort_fo2";
    # system "bedtools intersect -a $sort_fo2 -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -wo  |gzip > ../../../output/figure/whole_genome/3pca_1e_4/HIC/hotspot_cluster_${sample}_5000_inter_chr.bed.gz";
    # print "$sample\n";
# }







