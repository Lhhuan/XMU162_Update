#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

use List::MoreUtils ':all';

my $f1 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/HISTONE_MARK_AND_VARIANT/merge_pos_info_sample_narrow_peak.bed.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $marker = $f[2];
        my $DCid =$f[3];
        my $cell_line =$f[4];
        if($marker=~/HISTONE_MARK_AND_VARIANT/){
            # print "$DCid\n";
            $hash1{$DCid}=1;
        }
    }
    # my $v="$DCid\t$cell_line";
    # push@{$hash1{$marker}},$DCid;
}

while(<$I2>)
{
    chomp;
    unless(/^study/){
        my @f=split/\t/;
        my $id= $f[3];
        $hash2{$id}=1;
    }
}

foreach my $k(sort keys %hash1){
    unless(exists $hash2{$k}){
        print "$k\n";
    }
}
print "=======================";

foreach my $k(sort keys %hash2){
    unless(exists $hash1{$k}){
        print "$k\n";
    }
}