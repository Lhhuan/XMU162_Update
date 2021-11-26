
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/Cis_eQTL/06_merge_all_tissue_cis_sig_eQTL_hotspot_egene.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_sorted.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "./output/hotspot_without_egene.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;


while(<$I1>)
{
    chomp;
    unless(/^Chr/){
        my @f = split/\t/;
        my $k =join("\t",@f[0..2]);
        $hash1{$k}=1;
        # print "$k\n";
    }
}

while(<$I2>)
{
    chomp;
    my @f = split/\t/;
    my $k =join("\t",@f[0..2]);
    unless(exists $hash1{$k}){
        print $O1 "$_\n";
        # print "$_\n";
    }
}

close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ./output/hotspot_without_egene_sorted.bed.gz";




