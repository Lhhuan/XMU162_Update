#/home/huanhuan/project/eQTL_Catalogue/output/02_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted_chr${i}.bed.gz 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"进行 annotation得../output/edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz, 然后汇总调整格式得"../output/edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my %hash1;
my $sorted_input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
my $fo1 = "../output/edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

for (my $i=1;$i<23;$i++){
    print "$i\tstart\n";
    my $f1= "../output/edges_annotation/hotspot_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted_chr${i}.bed.gz";
    system "bedtools intersect -a $sorted_input_file -b /home/huanhuan/project/eQTL_Catalogue/output/02_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted_chr${i}.bed.gz -wa -wb |cut -f1-3,7-9|gzip > ${f1}" ;
    # print $O1 "hchr_hstart_hend\tegene\ttissue:p_value";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    while(<$I1>)
    {
        chomp;
        my @f=split/\t/;
        my $hotspot=join("_",@f[0..2]);
        my $ensg= $f[3];
        my $p_value= $f[4];
        my $tissue =$f[5];
        my $k = "$hotspot\t$ensg";
        my $v = "$tissue:$p_value";
        push @{$hash1{$k}},$v;
    }
    print "$i\tend\n";
}

foreach my $k(sort keys %hash1){
    my @vs=@{$hash1{$k}};
    @vs=uniq(@vs);
    my $v= join(";",@vs);
    print $O1 "$k\t$v\n";
}