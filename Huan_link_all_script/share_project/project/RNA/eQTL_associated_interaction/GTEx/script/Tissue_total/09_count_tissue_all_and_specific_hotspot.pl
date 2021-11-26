#得"../../output/Tissue_total/number_of_tissue_all_hotspot.txt.gz"，"../../output/Tissue_total/number_of_tissue_specific_hotspot.txt.gz"
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


 
my $fo1 = "../../output/Tissue_total/number_of_tissue_all_hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo3 = "../../output/Tissue_total/number_of_tissue_specific_hotspot.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;

print $O3 "Tissue\tNumber\n";
print $O1 "Tissue\tNumber\n";
 

foreach my $tissue(@tissues){
    my $tissue11 = $tissue; 
    $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
    my $all_hotspot = "../../output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $tissue_specific_hotspot = "../../output/Tissue_total/specific/pure/${tissue}_specific.bed.gz";
    my @arg1s = stat ($all_hotspot);
    my $all_hotspot_size = $arg1s[7];
    #------------------
    my @arg2s = stat ($tissue_specific_hotspot);
    my $hotspot_specific_size = $arg2s[7];
    #----------------------------------
    my $command_all_hotspot = "zless $all_hotspot | wc -l" ;
    my $all_hotspot_line_count = wc($command_all_hotspot);
    print $O1 "$tissue\t$all_hotspot_line_count\n";

    my $command_tissue_specific_hotspot = "zless $tissue_specific_hotspot  | wc -l" ;
    my $tissue_specific_hotspot_line_count = wc($command_tissue_specific_hotspot);
    print $O3 "$tissue\t$tissue_specific_hotspot_line_count\n";
    print "$tissue\n";
}


sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
