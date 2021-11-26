#计算"../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz" 中cutoff下每条染色体平均hotspot的数目和长度得"../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz"

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use POSIX;
use List::Util qw( min max );

my @tissues = ("Lung","Whole_Blood");
foreach my $tissue(@tissues){
    my (%hash1,%hash2,%hash3);
    my $f1 = "../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    #--------
    my $fo1 = "../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "Cutoff\tChr\tAverage_Length\tMax_Length\tMin_Length\tMax_Length/Min_Length\tMax_Length/Average_Length\thotspot_Number\n";
    #--------
    while(<$I1>){
        chomp;
        unless(/^Chr/){
            my @f = split/\t/;
            my $Chr = $f[0];
            my $Length = $f[1];
            my $Number = $f[2];
            my $Cutoff = $f[3];
            my $k= "$Cutoff\t$Chr";
            my $all_length = $Length *$Number;
            push @{$hash1{$k}},$Number;
            push @{$hash2{$k}},$all_length;
            push @{$hash3{$k}},$Length;

        }
    }

    foreach my $k(keys %hash1){
        my @nums = @{$hash1{$k}};
        my @all_lens = @{$hash2{$k}};
        my @lengths = @{$hash3{$k}};
        my $num  = sum @nums;
        my $len = sum @all_lens;
        my $max_length = max @lengths;
        my $min_length = min @lengths;
        my $max_min_times = $max_length/$min_length; 
        my $average_len = POSIX::ceil($len/$num); #向上取整
        my $max_mean_times = $max_length/$average_len;
        print $O1 "$k\t$average_len\t$max_length\t$min_length\t$max_min_times\t$max_mean_times\t$num\n";
    }
}


