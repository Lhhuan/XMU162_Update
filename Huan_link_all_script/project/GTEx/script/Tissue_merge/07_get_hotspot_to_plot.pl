#将每个组织进行显著eQTL-eGene的提取，得"$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz"，排序后和排序后的"../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot_sorted.txt.gz"进行 bedtools intersect,得
#"../../output/Tissue_total/gene/49_share_hotspot_${tissue}_gene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $fo1 = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region.bed";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $region_start = 144275000;
my $region_end = $region_start+ 10e+6;
print $O1 "chr1\t$region_start\t$region_end\n";
my $win_region="../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_100_win.bed.gz";
my  $hotspot ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz";
my $hotspot_in_region = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_hotspot.bed";
my $eqtl_gene = "../../output/Tissue_merge/Cis_eQTL/06_sig_eQTL_gene_TSS_sorted.bed.gz";
# my $eqtl_gene_no_header = "../output/Tissue_merge/Cis_eQTL/interaction_heatmap/06_merge_all_tissue_cis_sig_eQTL_hotspot_egene_no_header.txt.gz";
my $qtl_in_hotspot_region = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_hotspot_eqtl_egene.bed.gz";
my $qtl_in_win = "../../output/Tissue_merge/Cis_eQTL/interaction_heatmap/07_1MB_region_hotspot_eqtl_in_win.bed.gz";
# zless $eqtl_gene |awk 'NR >1 {print}'|gzip 
system "bedtools intersect -a $fo1 -b $hotspot -wo |cut -f4-6  >$hotspot_in_region";
system "bedtools intersect -a $hotspot_in_region -b $eqtl_gene -wo |gzip>$qtl_in_hotspot_region"; #hotspot snp tss

system "bedtools makewindows -b $fo1 -w 100 -i winnum | gzip >$win_region";
# system  "zless $qtl_in_hotspot_region |awk '{print $4 }'"
# system "bedtools intersect -a $qtl_in_hotspot_region -b $win_region -wo |gzip >$qtl_in_win";








# my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");

# my $out_dir ="/share/data0/GTEx/data/GTEx_v8_hg19_eQTL_egene";

# my $all_tissue_h = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz";
# # system "zless /share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz |sort -k1,1 -k2,2n |gzip > $all_tissue_h";

# my $fo1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/06_merge_all_tissue_cis_sig_eQTL_hotspot_egene.txt.gz";
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Chr\tStart\tEnd\tegene\thotspot\ttissue\n";


# foreach my $tissue(@tissues){
#     # print "$tissue\tstart\n";
#     my $tissue_eqtl_egene = "$out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz";
#     # close($O1);
#     system "bedtools intersect -a $all_tissue_h -b $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz -wo |cut -f1-3,7|sort -u > tmp.bed";
#     my $f1 = "tmp.bed";
#     open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
#     # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#     print "$tissue\tstart\n";
#     my %hash1;
#     while(<$I1>)
#     {
#         chomp;
#         my @f =split/\t/;
#         my $chr = $f[0];
#         my $start = $f[1];
#         my $end = $f[2];
#         my $egene = $f[3];
#         my $hotspot = join("_",@f[0..2]);
#         print $O1 "$_\t$hotspot\t$tissue\n";
#         # print "$_\n";
#     }
#     print "$tissue\tend\n";
# }



