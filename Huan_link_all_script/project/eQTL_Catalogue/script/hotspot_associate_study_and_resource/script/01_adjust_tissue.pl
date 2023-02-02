#need_study_for_hotspot_download_tabix_ftp_paths.tsv 中有些与gtex记录不相同，进行调整得../output/01_adjust_tissue_need_study_for_hotspot.tsv
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/eQTL_Catalogue/output/need_study_for_hotspot_download_tabix_ftp_paths.tsv" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "../output/01_adjust_tissue_need_study_for_hotspot.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";






while(<$I1>)
{
    chomp;
    if(/^study/){
        print $O1 "$_\n";
    }
    else{
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $tissue_ontology_id =$f[2];
        my $tissue_ontology_term =$f[3];
        my $tissue_label = $f[4];
        my $condition_label =$f[5];
        my $quant_method =$f[6];
        my $sample_size =$f[7];
        my $ftp_path =$f[8];
        my $refine_tissue_label =$f[9];
        print $O1 "$study\t$qtl_group\t$tissue_ontology_id\t$tissue_ontology_term\t$tissue_label\t$refine_tissue_label\t$condition_label\t$quant_method\t$sample_size\t$ftp_path\n";
    }
}
