#利用../../output/cancer_total/share/total/05_cancer_share_total.bed.gz得#../../output/cancer_total/11_count_cancer_tissue_share_hotspot.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $cancer = "Lung";
my $j = 18;


my $fo1 = "../../output/cancer_total/11_count_cancer_tissue_share_hotspot111.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Number\tcancer\ttissue\n";


my $f1 = "../../output/cancer_total/share/total/05_cancer_share_total.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $v =join("\t",@f[0..2]);
    my $cancer = $f[-2];
    my $tissue = $f[-1];
    my $k = "$cancer\t$tissue";
    push @{$hash1{$k}},$v;
}    

foreach my $k(sort keys %hash1 ){
    my @vs = @{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs; 
    my $count =@vs;
    print $O1 "$count\t$k\n";
}







sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
