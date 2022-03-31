#将所有tissue "${dir}/${tissue}${suffix}"合并得gene文件"../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz";
#合并位点并用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全得"../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");

my $fo2 = "../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_sQTL_gene.txt.gz";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O2, "| gzip >$fo2" or die $!;

my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL_hg19";
my $suffix = ".v8.sqtl_signifpairs.txt.gz";
# print $O1 "SNP_chr\tSNP_pos\tPvalue\n";
print $O2 "SNP_chr\tSNP_pos\tPvalue\tgene_id\n";
my (%hash1,%hash2,%hash3);

foreach my $tissue(@tissues){
    print "$tissue\n";
    # $hash1{$tissue}=1;
    my $f1 ="${dir}/${tissue}${suffix}";; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        while(<$I1>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $phenotype_id = $f[3];
            my $Pvalue =$f[-6];
            my $k1 = "$SNP_chr\t$SNP_pos\t$Pvalue";
            my @t=split/\:/,$phenotype_id;
            my $gene_id = $t[-1];
            my $k2 = "$k1\t$gene_id";
            unless(exists $hash2{$k2}){
                $hash2{$k2}=1;
                print $O2 "$k2\n";
            }
        }
    }
}
