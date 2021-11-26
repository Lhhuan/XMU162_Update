#将hotspot makewindows 1MB  ../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}_makewin_1MB.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;
# my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");

# my $command1= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  200000 | gzip > hg19.chrom1_22_sizes_sorted.20KB_windows.bed.gz";
# my $command1= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  10000000 | gzip > hg19.chrom1_22_sizes_sorted.10MB_windows.bed.gz";
# my $command1= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  1000 | gzip > hg19.chrom1_22_sizes_sorted.1KB_windows.bed.gz";
# my $command2= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  10000 | gzip > hg19.chrom1_22_sizes_sorted.10KB_windows.bed.gz";
# my $command2= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  5000 | gzip > hg19.chrom1_22_sizes_sorted.5KB_windows.bed.gz";
my $command2= "bedtools makewindows -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt -w  261 | gzip > hg19.chrom1_22_sizes_sorted.261BP_windows.bed.gz";
system "$command2";
my @tissues =  ("Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Prostate","Thyroid","Whole_Blood");
my $win = "261BP";
foreach my $tissue1(@tissues){
    my $tissue11 = $tissue1; 
    $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
    my $file1 = "../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz"; 
    my $command2 = "bedtools intersect -a $file1 -b hg19.chrom1_22_sizes_sorted.${win}_windows.bed.gz -wo |cut -f4-7 |gzip > ../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}_makewin_${win}.bed.gz";
    system $command2;

}
