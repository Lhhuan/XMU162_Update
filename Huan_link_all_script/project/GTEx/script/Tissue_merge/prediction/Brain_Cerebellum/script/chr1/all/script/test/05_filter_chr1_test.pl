#利用 "/home/huanhuan/reference/grch37_ensgID_position_from_ensembl104.txt" 和"../output/01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt"为"../output/04_transform_and_filter_predict_top0.1.txt" 提取hotspot 和 egene 的位置，得 "../output/05_all_predict_top0.1.bed.gz"，"../output/05_trans_predict_top0.1.bed.gz"，"../output/05_cis_predict_top0.1.bed.gz"，并得"../output/01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt" egene 位置得"../output/05_groundtruth_link_pos.bed.gz"

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/output/01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "/home/huanhuan/reference/grch37_ensgID_position_from_ensembl104.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $f3 = "/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/output/04_transform_and_filter_predict.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/05_trans_predict.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/05_cis_predict.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
my $fo3 = "../output/05_all_predict.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

my $fo4 = "../output/05_groundtruth_link_pos.bed.gz";
open my $O4, "| gzip >$fo4" or die $!;

my $header ="hotspot_id\tgene_id\tsmiliarity\th_chr\th_start\th_end\tgene\tgene_chr\tgene_start\tgene_end";
print $O1 "$header\n";
print $O2 "$header\n";
print $O3 "$header\n";
my (%hash1,%hash2,%hash3,%hash4);

while(<$I2>)
{
    chomp;
    if(/^ENSG/){
        my @f=split/\t/;
        my $ensg = $f[0];
        my $chr= $f[1];
        my $start =$f[2];
        $start = $start-1; #1-based
        my $end = $f[3];
        my $v = "$chr\t$start\t$end";
        if($chr =~/\b1\b/){
            $hash4{$ensg}=$v;
            # print "$ensg\n";
        }

    }
}


while(<$I1>)
{
    chomp;
    unless(/^Chr/){
        my @f = split/\t/;
        my $chr = $f[0];
        my $start =$f[1];
        my $end=$f[2];
        my $egene= $f[3];
        my $hotspot =$f[4];
        # print "$egene\t$hotspot\n";
        my $hotspotidx =$f[5];
        my $geneidx =$f[6];
        if($chr =~/\bchr1\b/){
            $hash1{$hotspotidx}=$hotspot;
            $hash2{$geneidx}=$egene;
            my $edge1 = "$hotspotidx\t$geneidx";
            my $edge2 = "$geneidx\t$hotspotidx";
            $hash3{$edge1}=1;
            $hash3{$edge2}=1;
            if(exists ($hash4{$egene})){
                my $gene_pos = $hash4{$egene};
                # print "$egene\n";
                print $O4 "$chr\t$start\t$end\t$egene\t$gene_pos\n"; #h_chr\th_start\th_end\tegene\tgene_chr\tgene_start\tgene_end\n
            }
        }
    }
}

while(<$I3>)
{
    chomp;
    unless(/^hotspot_id/){
        my @f=split/\t/;
        my $hotspot_id = $f[0];
        my $gene_id= $f[1];
        my $smiliarity =$f[2];
        my $k = "$hotspot_id\t$gene_id";
        unless(exists $hash3{$k}){   
            if(exists $hash1{$hotspot_id}){
                my $hotspot = $hash1{$hotspot_id};
                my @h=split/_/,$hotspot;
                my $h_chr = $h[0];
                my $h_start=$h[1];
                my $h_end =$h[2];
                my $gene = $hash2{$gene_id};
                # $gene =~ s/.*//g;
                if(exists $hash4{$gene}){
                    my $gene_pos =$hash4{$gene};
                    my @t=split/\t/,$gene_pos;
                    my $gene_chr = $t[0];
                    my $gene_start = $t[1];
                    my $gene_end =$t[2];
                    $gene_chr="chr${gene_chr}";
                    print $O3 "$_\t$h_chr\t$h_start\t$h_end\t$gene\t$gene_chr\t$gene_start\t$gene_end\n";
                    if($gene_chr ne $h_chr){
                        print $O1 "$_\t$h_chr\t$h_start\t$h_end\t$gene\t$gene_chr\t$gene_start\t$gene_end\n";
                    }
                    else{
                        print $O2 "$_\t$h_chr\t$h_start\t$h_end\t$gene\t$gene_chr\t$gene_start\t$gene_end\n";
                    }
                }
            }
        }
    }
}


# close($I1);
# close($O1);