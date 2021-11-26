 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");
# my @exists_tissue = ("Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");
# my @un_exists_tissue =("Kidney_Medulla");

my %hash1;
foreach my $tissue(@tissues){
    $hash1{$tissue}=1;
}

my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
opendir (DIR, $dir) or die "can't open the directory!";
my  @dir = readdir DIR;
foreach my $file(@dir){
    if($file =~$suffix){

        $file =~ s/$suffix//g;
        unless(exists $hash1{$file}){
            print "$file\n";
        }
    }
}
