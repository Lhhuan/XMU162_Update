#将每个组织进行显著eQTL-snp的chr1进行提取，得"./tmp_output/07_chr_muliti_tissue_eqtl.bed.gz， muti-tissue hotspot的chr1进行提取得 ./tmp_output/07_chr_hotspot.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");

# my $out_dir ="/share/data0/GTEx/data/GTEx_v8_hg19_eQTL_egene";

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_sorted.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "./tmp_output/07_chr_hotspot.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "./tmp_output/07_chr_muliti_tissue_eqtl.bed.gz";
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O2, "| gzip >$fo2" or die $!;
# print "$tissue\tstart\n";
my %hash1;
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $chr = $f[0];
    my $start = $f[1];
    my $end = $f[2];
    if($chr=~ /\bchr1\b/){
    # if($chr=~ /chr1/){
        print $O1 "$_\n";
    }
}

my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";

foreach my $tissue(@tissues){
    print "$tissue\n";
    my $f2 ="${dir}/${tissue}${suffix}";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I2>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $variant_id = $f[0];
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $Pvalue =$f[-6];
            my $start = $SNP_pos;
            my $end = $SNP_pos +1;
            my @t = split/\_/,$variant_id;
            my $ref = $t[2];
            my $alt = $t[3];
            # print "$gene_id1\t$gene_id\n";
            if($Pvalue < 5E-8 && $SNP_chr =~/\b1\b/){
                my $output = "chr$SNP_chr\t$start\t$end\t$ref\t$alt\t$tissue";
                unless(exists $hash1{$output}){
                    $hash1{$output}=1;
                    print $O2 "$output\n";
                }
            }
        }
    }
}

close($O2);

# system "zless $fo2 |sort -k1,1 -k2,2n |gzip > ./tmp_output/01_chr_hotspot_sorted.bed.gz";
