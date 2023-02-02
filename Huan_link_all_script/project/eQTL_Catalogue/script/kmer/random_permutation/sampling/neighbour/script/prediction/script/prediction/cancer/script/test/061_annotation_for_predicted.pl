
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $sorted_input_file = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/warm_hotspot_region_win4518_large_than6.bed.gz";
my $input_file_base_name = basename($sorted_input_file);
my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/ca_tf_TCGA_sample.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^marker/){
        my @f=split/\t/;
        my $marker =$f[0];
        my $TCGA =$f[1];
        my $k ="$marker\t$TCGA";
        $hash1{$k}=1;
    }
}

foreach my $k (sort keys %hash1){
    print "$k\n";
    my @f =split/\t/,$k;
    my $marker =$f[0];
    my $TCGA =$f[1];
    my $marker1 = $marker;
    $marker1 =~s/Human_CHROMATIN_Accessibility/CHROMATIN_Accessibility/g;
    $marker1 =~s/Human_FACTOR/TFBS/g;
    my $annotation_file="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/${marker}/${TCGA}/merge_mean_signalvalue.bedgraph.gz";
    my $out_dir ="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/cancer/output/${TCGA}";
    my $out_file = "$out_dir/${marker1}_${input_file_base_name}";
    unless(-e $out_dir ){
        system "mkdir -p $out_dir ";
    }
    system "bedtools intersect -a $sorted_input_file -b $annotation_file -wo |gzip >$out_file";
}

#===============

my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/TCGA_mark.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my %hash2;
while(<$I2>)
{
    chomp;
    my @f=split/\t/;
    my $TCGA =$f[0];
    my $marker =$f[1];
    my $k ="$marker\t$TCGA";
    $hash2{$k}=1;
}


foreach my $k (sort keys %hash2){
    print "$k\n";
    my @f =split/\t/,$k;
    my $marker =$f[0];
    my $TCGA =$f[1];
    my $annotation_file="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/HISTONE_MARK_AND_VARIANT/${TCGA}/${marker}/merge_mean_signalvalue.bedgraph.gz";
    my $out_dir ="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/cancer/output/${TCGA}";
    my $out_file = "$out_dir/${marker}_${input_file_base_name}";
    unless(-e $out_dir ){
        system "mkdir -p $out_dir ";
    }
    system "bedtools intersect -a $sorted_input_file -b $annotation_file -wo |gzip >$out_file";
}