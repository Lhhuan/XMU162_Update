#将每个组织进行显著eQTL-eGene的提取，得"$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz"，排序后和排序后的"../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot_sorted.txt.gz"进行 bedtools intersect,得
#"../../output/Tissue_total/gene/49_share_hotspot_${tissue}_gene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use Env qw(PATH);
my $region = "half_chr1";
my $win =100000;
my $fo1 = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}.bed";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $region_start = 125000000;
my $region_end = 242249719;
print $O1 "chr1\t$region_start\t$region_end\n";
my $win_region="../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}_win${win}.bed.gz";
# my  $hotspot ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz";
my  $hotspot ="/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed";
my $hotspot_in_region = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}_hotspot.bed";
my $eqtl_gene = "../../output/Tissue_merge/Cis_eQTL/06_sig_eQTL_gene_TSS_sorted.bed.gz";
# my $eqtl_gene_no_header = "../output/Tissue_merge/Cis_eQTL/interaction_heatmap/06_merge_all_tissue_cis_sig_eQTL_hotspot_egene_no_header.txt.gz";
my $qtl_in_hotspot_region = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}_hotspot_eqtl_egene.bed.gz";
my $qtl_in_r =  "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_${region}_eQTL.bed.gz";
my $tss_in_r =  "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_in_region_${region}_eQTL_TSS.bed.gz";
my $qtl_in_win = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}_hotspot_eqtl_in_win${win}.bed.gz";
my $tss_in_win = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_region_${region}_hotspot_eqtl_TSS_in_win${win}.bed.gz";
# zless $eqtl_gene |awk 'NR >1 {print}'|gzip 
# system "bedtools intersect -a $fo1 -b $hotspot -wo |cut -f4-6  >$hotspot_in_region";
# system "bedtools intersect -a $hotspot_in_region -b $eqtl_gene -wo |gzip>$qtl_in_hotspot_region"; #hotspot snp tss

# system "bedtools makewindows -b $fo1 -w $win -i winnum | gzip >$win_region";
# system  "zless $qtl_in_hotspot_region |awk -v OFS='\t' '{print $4,$5,$6,$7}' |sort -u |sort -k1,1 -k2,2n | gzip > $qtl_in_r";
# system  "zless $qtl_in_hotspot_region |awk -v OFS='\t' '{print $8,$9,$10,$7}' |sort -u |sort -k1,1 -k2,2n | gzip > $tss_in_r";


# system "bedtools intersect -a $qtl_in_r -b $win_region -wa -wb |cut -f1-4,8|gzip > $qtl_in_win";
# system "bedtools intersect -a $tss_in_r -b $win_region -wa -wb |cut -f1-4,8|gzip > $tss_in_win";
$ENV{'fo1'} = $fo1 ;
$ENV{'hotspot'} = $hotspot ;
$ENV{'eqtl_gene'} = $eqtl_gene ;
$ENV{'hotspot_in_region'}  = $hotspot_in_region;
$ENV{'qtl_in_hotspot_region'} = $qtl_in_hotspot_region ;
$ENV{'win'} = $win ;
$ENV{'qtl_in_r'} = $qtl_in_r ;
$ENV{'win_region'} = $win_region ;
$ENV{'tss_in_r'} = $tss_in_r ;
$ENV{'qtl_in_win'} = $qtl_in_win ;
$ENV{'tss_in_win'} = $tss_in_win ;

 #设置环境变量
system "bash bedtools_process.sh"




