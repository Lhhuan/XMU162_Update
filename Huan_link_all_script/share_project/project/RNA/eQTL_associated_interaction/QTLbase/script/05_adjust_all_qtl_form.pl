#调整../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz 的格式，得
# ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

#------------------------------------

my $fo1 = "../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz ";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
# my $fo2 = "../output/emplambda_in_chr_max_position.txt";
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

# print $O1 "chr\tmin\tmax\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $t =$f[1]; 
        my $pos = $t; #snp pos
        # $pos =~ s/"//g;
        my $start = $pos-1;
        my $end = $pos;
        my $chr =$f[2]; #snp chr
        my $xQTL = $f[3];
        print $O1 "$chr\t$start\t$end\t$xQTL\t$emplambda\n";
    }
}
