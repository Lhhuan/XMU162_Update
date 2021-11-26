#根据"../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz"和"/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes" split 染色体并random取hotspot的均长，得../


#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use POSIX;
use List::Util qw( min max );

my $f1 = "/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my %hash1;
while(<$I1>){
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];
    my $size = $f[1];
    $hash1{$chr}=$size;
}


my @tissues = ("Lung");
# my @tissues = ("Lung","Whole_Blood");
my @cutoffs;
push @cutoffs,0.01;
for (my $i=0.05;$i<0.7;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}

# foreach my $cutoff()

foreach my $tissue(@tissues){
    my (%hash2,%hash3);

    foreach my $Cutoff(@cutoffs){
    my $f2 = "../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
        # print "$Cutoff\n";
        while(<$I2>){
            chomp;
            unless(/^Cutoff/){
                my @f = split/\t/;
                my $cutoff =$f[0];
                my $Chr = $f[1];
                my $Average_Length = $f[2];
                my $hotspot_Number = $f[-1];
                # print "$cutoff\n";
                if( $Chr =~/\bchr1\b/ ||$Chr =~/\bchr22\b/ ){
                    if(abs($cutoff -$Cutoff)<0.0000001){
                    # if($cutoff ==$Cutoff){
                        print "$cutoff\t$Cutoff\t$Chr\t$Average_Length\n";
                    }
                }
            }
        }
    }
}
