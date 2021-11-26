
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


my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");


my $f1 = "../../output/Tissue_total/share/all/07_all_tissue_intersect.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件  

my $fo3 = "../../output/Tissue_total/specific/pure/tissue_specific.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

print $O3 "h_chr\th_start\th_end\ttissue\n";

#----------------





my %hash2;
while(<$I1>)
{
    chomp;
    unless(/^tissue1_chr/){
        my @f= split/\t/;
        my $t1_chr =$f[0];
        my $t1_start =$f[1];
        my $t1_end = $f[2];
        my $t2_chr = $f[3];
        my $t2_start =$f[4];
        my $t2_end = $f[5];
        my $k2 = "$t1_chr\t$t1_start\t$t1_end";
        my $k3 = "$t2_chr\t$t2_start\t$t2_end";
        $hash2{$k2}=1;
        $hash2{$k3}=1;
    }
}    

foreach my $tissue(@tissues){
    my $tissue11 = $tissue; 
    $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
    my $file1 = "../../output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $f3 = $file1; #cancer
    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件  

    my $fo1 = "../../output/Tissue_total/specific/pure/${tissue}_specific.bed.gz";
    open my $O1, "| gzip >$fo1" or die $!;


    while(<$I3>)
    {
        chomp;
        my @f= split/\t/;
        my $chr =$f[0];
        my $start =$f[1];
        my $end = $f[2];
        my $k = "$chr\t$start\t$end";
        unless(exists $hash2{$k}){
            print $O1 "$_\n";
            print $O3 "$_\t$tissue\n";
        }
    } 
    print "$tissue\n";
}
