#对 "../../output/Tissue_total/share/total/05_absolute_tissue_intersect.bed.gz" 进行 tissue合并得../../output/Tissue_total/share/total/06_all_tissue_share_hotspot_total_contain.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;

my $f1 = "../../output/Tissue_total/share/total/05_absolute_tissue_intersect.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


my $fo1 = "../../output/Tissue_total/share/total/06_all_tissue_share_hotspot_total_contain.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;


print $O1 "chr\tstart\tend\tshare_tissues\n";

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^tissue1_chr/){
        my @f= split/\t/;
        my $tissue1_chr = $f[0];
        my $tissue1_start =$f[1];
        my $tissue1_end =$f[2];
        my $tissue1 =$f[-2];
        my $tissue2 =$f[-1];
        my $k = "$tissue1_chr\t$tissue1_start\t$tissue1_end";
        push @{$hash1{$k}},$tissue1;
        push @{$hash1{$k}},$tissue2;
        # print $O1 "$_\t$tissue1\t$tissue2\n";
    }
}   

foreach my $k(sort keys %hash1){
    my @vs = @{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs; 
    my @sorted_vs = sort {$b cmp $a} @vs;
    my $share_tissues = join (",",@sorted_vs);
    print $O1 "$k\t$share_tissues\n";
}
