#对 "../../output/cancer_total/share/total/17_absolute_cancer_cancer_intersect.bed.gz" 进行 tissue合并得"../../output/cancer_total/share/total/18_all_tissue_cancer_hotspot_total_contain.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;

my $f1 = "../../output/cancer_total/share/total/17_absolute_cancer_cancer_intersect.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


my $fo1 = "../../output/cancer_total/share/total/18_all_tissue_cancer_hotspot_total_contain.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;


print $O1 "chr\tstart\tend\tshare_cancers\n";

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^cancer1_chr/){
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
