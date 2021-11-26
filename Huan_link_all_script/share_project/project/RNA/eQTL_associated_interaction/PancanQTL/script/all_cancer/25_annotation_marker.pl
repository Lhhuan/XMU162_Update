#annotation marker ,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/HISTONE_MARK_AND_VARIANT/${cancer}/${marker}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;
my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");
# my @cancers =  ("KIRP"); 
my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/TCGA_mark.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  
my %hash1;
while(<$I2>)
{
    chomp;
    my @f = split/\t/;
    my $cancer = $f[0];
    my $marker = $f[1];
    push @{$hash1{$cancer}},$marker;
}



my $anno_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell";
foreach my $cancer(@cancers){
    # my $f1 = "../../output/${cancer}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${cancer}.txt.gz";
    my $f1 = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $out_dir = "../../output/${cancer}/Cis_eQTL/anno";
    # system "rm -r $out_dir";
    if(-e $out_dir){
    print "${out_dir}\texist\n";
    }
    else{
        system "mkdir -p $out_dir";
    }
    my $input_file_base_name = basename($f1);
    my $tfbs_file = "$anno_dir/Human_FACTOR/${cancer}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";
    my $tfbs = "bedtools intersect -a $f1 -b $tfbs_file -wo |gzip >${out_dir}/TFBS_${input_file_base_name}";
    if (-e $tfbs_file){ #判断文件存在
        system "$tfbs";
    }
    my $ca_file = "$anno_dir/Human_CHROMATIN_Accessibility/${cancer}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";
    my $ca = "bedtools intersect -a $f1 -b $ca_file  -wo |gzip >${out_dir}/CHROMATIN_Accessibility_${input_file_base_name}";
    if(-e $ca_file){
        system "$ca";
    }
    #-----------------------
    if(exists $hash1{$cancer}){
        my @markers = @{$hash1{$cancer}};
        foreach my $marker(@markers){
            # "./cancer_cell/${out_dir}/${TCGA}/${factor}";
            my $marker_file = "$anno_dir/HISTONE_MARK_AND_VARIANT/${cancer}/${marker}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";
            my $marker_commmand = "bedtools intersect -a $f1 -b $marker_file -wo |gzip >${out_dir}/${marker}_${input_file_base_name}";
            system "$marker_commmand ";
        }
    }
}





