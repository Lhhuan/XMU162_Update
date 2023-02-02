#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
my $f1 = "./output/02_tissue_level_Marker_source.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info_tissue_refine.tsv";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^tissue_ontology_id/){
        my @f = split/\t/;
        my $refine_tissue_label2 =$f[1];
        $hash1{$refine_tissue_label2}=1;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $refine_tissue_label2 = $f[6];
        my $ftp_path =$f[10];
        my @t= split/\//,$ftp_path;
        my $filename=$t[-1];
        # print "$study\t$refine_tissue_label2\t$filename\n";
        my $v = "$study/$filename";
        if(exists $hash1{$refine_tissue_label2}){
            push @{$hash2{$refine_tissue_label2}},$v;
        }
    }
}

my $eqtl_dir="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data";
foreach my $tissue(sort keys %hash2){
    # my $pid = $pm->start and next; #开始多线程
    print "$tissue\tstart1\n";
    my @vs=@{$hash2{$tissue}};
    my $fo1 = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/08_eqtl.bed.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    foreach my $v(@vs){
        print "$tissue\t$v\n";
        my $f3 = "$eqtl_dir/$v";
        open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
        while(<$I3>)
        {
            chomp;
            unless(/^molecular_trait_id/){
                my @f = split/\t/;
                my $chr = $f[1];
                my $pos= $f[2];
                my $ref= $f[3];
                my $alt= $f[4];
                # my $variant =$f[5];
                my $pvalue =$f[8];
                my $gene_id = $f[-3];
                my $chrpos ="$chr\t$pos";
                my $variant = "$chrpos\t$ref\t$alt";
                my $start = $pos-1+1;
                my $end =$pos+1;
                
                my $bed="chr${chr}\t$start\t$end\t$pvalue\t$gene_id";
                print $O1 "$bed\n";
            }
        }
    }
    close($O1);
    my $fo2 = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/08_eqtl_sorted.bed.gz";
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip >$fo2";
    print "$tissue\tfinish1\n";
    # $pm->finish;  #多线程结束
}

print "output\tfinish\n";
#===================多线程输出，但线程排序
# foreach my $tissue(sort keys %hash2){
#     print "$tissue\tstart\n";
#     my $fo1 = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/08_eqtl.bed.gz";
#     my $fo2 = "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/08_eqtl_sorted.bed.gz";
#     system "zless $fo1 |sort -k1,1 -k2,2n |gzip >$fo2";
#     print "$tissue\tfinish\n";
# }