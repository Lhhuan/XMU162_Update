#对$input_dir/${mark}_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz 进行jaccard index 进行计算,得的marker进行计算,得"$out_dir/${cutoff2}_jaccard_index_histone_marker.txt.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $cutoff2 = "0_0.176";

my $input_dir ="/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/annotation/extend/filter_3103833";
my $input_file_base_name= "Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz";

my @markers = ("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","CHROMATIN_Accessibility","TFBS","CTCF");

my $fo1 = "../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotation_fraction.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Marker\tchr\tstart\tend\toverlap_fraction\n";

foreach my $mark (@markers){
    my %hash1;
    my $f1 = "$input_dir/${mark}_${input_file_base_name}";
    # print "$f2\n";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $hotspot_chr = $f[0];
        my $hotspot_start = $f[1];
        my $hotspot_end = $f[2];
        my $overlap_bp = $f[-1];
        my $overlap_fraction= $overlap_bp/($hotspot_end - $hotspot_start);
        my $k = "$hotspot_chr\t$hotspot_start\t$hotspot_end";
        $hash1{$k}=1;
        print $O1 "$mark\t$k\t$overlap_fraction\n";
    }
    my $f2= "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        my $hotspot_chr = $f[0];
        my $hotspot_start = $f[1];
        my $hotspot_end = $f[2];
        my $k ="$hotspot_chr\t$hotspot_start\t$hotspot_end";
        unless(exists $hash1{$k}){
            print $O1 "$mark\t$_\t0\n";
        }
        
    }
}

