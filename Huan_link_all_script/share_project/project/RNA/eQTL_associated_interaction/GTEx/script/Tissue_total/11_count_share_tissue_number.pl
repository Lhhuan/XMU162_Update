#统计 #"../../output/Tissue_total/share/total/06_all_tissue_share_hotspot_total_contain.bed.gz"得"../../output/Tissue_total/11_share_tissue_number_count.txt.gz"
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



my $f1 = "../../output/Tissue_total/share/total/06_all_tissue_share_hotspot_total_contain.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件  

my $fo1 = "../../output/Tissue_total/11_share_tissue_number_count.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "number_of_share_tissue\tnumber_of_hotspots\n";

#----------------





my %hash2;
while(<$I1>)
{
    chomp;
    # unless(/^chr$/){
        # print "$_\n";
    my @f= split/\t/;
    my $chr =$f[0];
    my $start =$f[1];
    my $end = $f[2];
    unless($start =~/start/){
        my $share_tissues = $f[3];
        my @t = split/,/,$share_tissues;
        my $length = @t;
        my $k = join("\t",@f[0..2]);
        # my $hash2{$length}
        push @{$hash2{$length}},$k;
    }
}

foreach my $k(sort keys %hash2){
    my @vs =@{$hash2{$k}};
    my %hash3;
    @vs = grep { ++$hash3{$_} < 2 } @vs; 
    my $number =@vs;
    print $O1 "$k\t$number\n";
    if($k == 49){
        my $v =join("\n",@vs);
        print "$v\t$k\n";
    }
    # print $O1 "$k\t$number";
}

