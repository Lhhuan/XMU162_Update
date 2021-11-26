#download data 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @array = ("786-O", "A375", "A549", "Adipocyte", "AML_blast", "Astrocyte", "BE2C", "BJ", "Bronchia_epithelial", "B_cell_blood", "C4-2", "Caco-2", "Calu-3", "CCRF-CEM", "CD133+", 
"CD14+", "CD14+_monocyte", "CD19+", "CD20+", "CD34+", "CD36+", "CD4+", "CD8+", "Cerebellum", "CMK", "CNCC", "Colo320", "Colo829", "CUTLL1", "CyT49", "Denditric_cell", "DLD1", "DOHH2", 
"ECC-1", "Endometrial_stromal_cell", "endometrioid_adenocarcinoma", "ESC", "ESC_neuron", "ESC_NPC", "Esophagus", "EWS502", "Fetal_brain", "Fetal_heart", "Fetal_kidney", "Fetal_lung", 
"Fetal_muscle_leg", "Fetal_placenta", "Fetal_small_intestine", "Fetal_spinal_cord", "Fetal_stomach", "Fetal_thymus", "Fibroblast_foreskin", "Foreskin_keratinocyte", "FT246", "FT33", 
"GC_B_cell", "Gliobla", "GM10847", "GM12878", "GM12891", "GM12892", "GM18486", "GM18505", "GM18507", "GM18508", "GM18516", "GM18522", "GM18526", "GM18951", "GM19099", "GM19141", "GM19193", 
"GM19238", "GM19239", "GM19240", "H1", "H128", "H2171", "H54", "H9", "HACAT", "HCASMC", "HCC1954", "HCT116", "Heart", "HEK293", "HEK293T", "HeLa-S3", "Hela", "Hepatocyte", "HepG2", "HFF", 
"HL-60", "hMADS-3", "HMEC", "hNCC", "HSC", "HSMM", "HSMMtube", "HT1080", "HT29", "HuCCT1", "HUES64", "HUVEC", "IMR90", "iPSC", "Jurkat", "K562", "Kasumi-1", "KATO3", "KB", "KELLY", 
"Keratinocyte", "Kidney", "Kidney_cortex", "Left_ventricle", "LHCN-M2", "Liver", "LNCaP-1F5", "LNCaP-abl", "LNCaP", "LoVo", "LP-1", "LS174T", "Lung", "LY1", "Macrophage", "MCF-7", "MCF10A", 
"MDA-MB-231", "ME-1", "Melanocyte", "melanoma", "Mesendoderm", "MM1S", "Monocyte", "MS1", "MSC_BM", "Myotube", "Namalwa", "NB4", "NCCIT", "NGP", "NH-A", "NHBE", "NHDF", "NHEK", "NHLF", "NKC", 
"NT2-D1", "OCI-LY1", "OCI-Ly7", "Osteobl", "Osteoblast", "Ovary", "P493-6", "PANC-1", "Pancreas", "Pancreatic_islet", "PBMC", "PC3", "Plasma_cell_myeloma", "PrEC", "Raji", "Ramos", "REH", 
"Retina", "RPE", "RPTEC", "RS4-11", "SEM", "SGBS_adipocyte", "SH-SY5Y", "SK-MEL-5", "SK-N-MC", "SK-N-SH", "SK-N-SH_RA", "Skeletal_muscle", "SkMC", "Small_intestine", "Sperm", "Spleen", 
"T47D-MTVL", "T47D", "T98G", "TC-797", "Th1", "Th2", "Thymus", "Treg_cell", "Trophoblast", "U2OS", "U87", "Urothelial_cell", "VCaP", "ZR75-1", "ZR75-30");



for  my $tissue (@array){ 
my $link = "http://www.enhanceratlas.org/data/download/enhancer/hs/${tissue}.bed";
my $commod = "wget ${link}\n";
system "$commod";
}