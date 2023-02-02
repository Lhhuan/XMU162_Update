#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::MoreUtils ':all';
use File::Path;

my %hash1;
# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
my $f1 = "./cancer_cell/histone_TCGA_mark.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    unless(/^TCGA/){
        my $TCGA = $f[0];
        my $factor =$f[1];
        my $cell_line = $f[2];
        my $Dcid =$f[3];
        my $k = "$TCGA\t$factor\t$cell_line";
        push@{$hash1{$k}},$Dcid;
    }
}

my $fo1 = "./cancer_cell/092_final_histone_TCGA_mark.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "TCGA\tfactor\tcell_line\tcell_line_refine\tDcid\n";
foreach my $k(sort keys %hash1){
    print "$k\tstart\n";
    my @t=split/\t/,$k;
    my $TCGA = $t[0];
    my $marker =$t[1];
    my $cell_line =$t[2];
    my $dir = "HISTONE_MARK_AND_VARIANT/human_hm";
    my $cell_line1 = $cell_line;
    $cell_line1 =~s/\s+/_/g;
    $cell_line1 =~s/;/_/g;
    $cell_line1 =~s/\./_/g;
    $cell_line1 =~s/\'//g;
    my $outdir = "./cancer_cell/cell_type_level/$TCGA/$cell_line1";
    if (-e $outdir){
        print "$outdir exist\n";
    }
    else{
        mkpath($outdir);
    }
    #========================
    my $fo2 = "${outdir}/${marker}_merge_pos_info_narrow_peak.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    my $fo3 = "${outdir}/${marker}_merge_pos_info_narrow_peak_signalValue.bed.gz";
    open my $O3, "| gzip >$fo3" or die $!; 
    my @vs = @{$hash1{$k}};
    my @uniq_vs = uniq(@vs);
    my %hash2;
    foreach my $DCid(@uniq_vs){
        my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/$dir/${DCid}_sort_peaks.narrowPeak.bed.gz"; 
        if(-e $f2){
            my %hash3;
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            while(<$I2>)
            {
                chomp;
                my @f =split/\t/;
                my $chr = $f[0];
                my $start =$f[1];
                my $end = $f[2];
                my $signalValue =$f[6];
                print $O2 "$chr\t$start\t$end\n";
                print $O3 "$chr\t$start\t$end\t$signalValue\n";
                $hash2{$f2}=1;
                $hash3{$f2}=1
            }
            my $L3=keys %hash3;
            if($L3 >0){ #判断 $f2 是否空
                print $O1 "$TCGA\t$marker\t$cell_line\t$cell_line1\t$DCid\n";
            }
            # print "$L3\n";
        }
    }
    close($O2);
    close($O3);
    my $L2 =keys %hash2;
    if($L2 >0){ #判断$fo3 是否空
        my $fo4="${outdir}/${marker}_merge_pos_info_narrow_peak_sorted.bed.gz";
        my $fo5="${outdir}/${marker}_merge_pos_info_narrow_peak_signalValue_sorted.bed.gz";
        my $fo6="${outdir}/${marker}_merge_pos_info_narrow_peak_sorted_merge.bed.gz";
        system "zless $fo2 |sort -k1,1 -k2,2n |gzip >$fo4";
        system "zless $fo3 |sort -k1,1 -k2,2n |gzip >$fo5";
        system "bedtools merge -i $fo4 |gzip >$fo6";
        print "$k\tfinish\n";
    }
    else{
        system "rm $fo2";
        system "rm $fo3";
    }
}
