 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @tissues =  ("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate");
my @exists_tissue = ("Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood");
my @un_exists_tissue =("Kidney_Medulla");




my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
my $pm = Parallel::ForkManager->new(10);
foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    print "$tissue\tstart\n";
    my $f1 ="${dir}/${tissue}${suffix}";; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $f2 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    #---------------------------mkdir
    # $tissue =~s/-/_/g;
    my $output_dir = "../../output/${tissue}";
    unless(-e $output_dir){
        system "mkdir -p $output_dir";
    }
    #---------------
    my $fo1 = "$output_dir/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "SNP_chr\tSNP_pos\tPvalue\n";


    my %hash1;
    while(<$I1>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $Pvalue =$f[-6];
            my $k = "$SNP_chr\t$SNP_pos";
            $hash1{$k}=1;
            print $O1 "$k\t$Pvalue\n";
        }
    }


    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
            }
        }
    }
    print "$tissue\tend\n";
    $pm->finish;  #多线程结束
}

