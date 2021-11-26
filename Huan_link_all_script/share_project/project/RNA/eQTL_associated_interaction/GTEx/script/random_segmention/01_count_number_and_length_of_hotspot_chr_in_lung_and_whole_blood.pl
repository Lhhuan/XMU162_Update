#统计/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.g中每条染色体的hotspot的长度及对应的数目，得../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_lung_and_whole_blood.txt.gz  ../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo2 = "../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_lung_and_whole_blood.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;

print $O2 "Tissue\tChr\tLength\tNumber\tCutoff\n";

my %hash1;
my @cutoffs=();
push @cutoffs,0.01;

for (my $i=0.05;$i<0.7;$i=$i+0.05){
    push @cutoffs,$i;
}

my @tissues = ("Lung");

foreach my $tissue (@tissues){
    my $fo1 = "../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "Chr\tLength\tNumber\tCutoff\n";
    foreach my $cutoff(@cutoffs){
        my(%hash2,%hash3,%hash4);
        my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz";
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

        while(<$I1>){
            chomp;
            my @f = split/\t/;
            my $chr = $f[0];
            my $start = $f[1];
            my $end = $f[2];
            my $length =$end - $start;  #0-based,实际end是此时end-1,因为后面还有取片段+1，故此次不再进行-1操作
            my $k = "$chr\t$length";
            push @{$hash2{$k}},$start;
        }

        foreach my $k(sort keys %hash2){
            my @poss = @{$hash2{$k}};
            my %hash; 
            @poss = grep { ++$hash{$_} <2 } @poss;
            my $number = @poss;
            print $O1 "$k\t$number\t$cutoff\n";
            print $O2 "$tissue\t$k\t$number\t$cutoff\n";
            # foreach my $p(@poss){
            #     print "$tissue\t$k\t$p\n";
            # }
        }
    }
}



my @tissue2s = ("Whole_Blood");
foreach my $tissue (@tissue2s){
    my $fo1 = "../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "Chr\tLength\tNumber\tCutoff\n";
    foreach my $cutoff(@cutoffs){
        my(%hash2,%hash3,%hash4);
        my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz";
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

        while(<$I1>){
            chomp;
            my @f = split/\t/;
            my $chr = $f[0];
            my $start = $f[1];
            my $end = $f[2];
            my $length =$end - $start;  #0-based,实际end是此时end-1
            my $k = "$chr\t$length";
            push @{$hash2{$k}},$start;
        }

        foreach my $k(sort keys %hash2){
            my @poss = @{$hash2{$k}};
            my %hash; 
            @poss = grep { ++$hash{$_} <2 } @poss;
            my $number = @poss;
            print $O1 "$k\t$number\t$cutoff\n";
            print $O2 "$tissue\t$k\t$number\t$cutoff\n";
        }
    }
}