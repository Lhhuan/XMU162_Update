#将每个组织进行显著eQTL-eGene的提取，得"$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz"，排序后和排序后的"../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot_sorted.txt.gz"进行 bedtools intersect,得
#"../../output/Tissue_total/gene/49_share_hotspot_${tissue}_gene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");

my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";

my $out_dir ="/share/data0/GTEx/data/GTEx_v8_hg19_eQTL_egene";
my $all_tissue_h = "../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot_sorted.txt.gz";
system "zless ../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot.txt.gz |sort -k1,1 -k2,2n |gzip > $all_tissue_h";
my $anno_out_dir = "../../output/Tissue_total/gene";
my $pm = Parallel::ForkManager->new(10);
foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    print "$tissue\tstart\n";
    my $f1 ="${dir}/${tissue}${suffix}";; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $fo1 = "$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    # print $O1 "SNP_chr\tSNP_pos\tPvalue\n";


    my %hash1;
    while(<$I1>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $gene_id1 = $f[3];
            my $Pvalue =$f[-6];
            my $start = $SNP_pos;
            my $end = $SNP_pos +1;
            my $gene_id = $gene_id1;
            $gene_id =~ s/\..*+//g;
            # print "$gene_id1\t$gene_id\n";
            if($Pvalue < 5E-8){
                my $output = "chr$SNP_chr\t$start\t$end\t$gene_id";
                unless(exists $hash1{$output}){
                    $hash1{$output}=1;
                    print $O1 "$output\n";
                }
            }
        }
    }
    close($O1);
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip > $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz";
    system "bedtools intersect -a $all_tissue_h -b $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz -wo |cut -f1-3,7|sort -u |gzip > $anno_out_dir/49_share_hotspot_${tissue}_gene.txt.gz";
    # system "bedtools intersect -a $all_tissue_h -b $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz -wo |sort -u |gzip > $anno_out_dir/49_share_hotspot_${tissue}_gene.txt.gz" ;
    print "$tissue\tend\n";
    $pm->finish;  #多线程结束
}



