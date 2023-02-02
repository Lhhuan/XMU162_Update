#判断"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" 中来自不同study 的相同tissue id是否完全相同
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::MoreUtils ':all';


my $f1 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# my $fo1 = "./output/01_adjust_tissue_need_study_for_hotspot.tsv";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";




my %hash1;

while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $tissue_ontology_id =$f[2];
        my $tissue_ontology_term =$f[3];
        my $tissue_label = $f[4];
        my $refine_tissue_label =$f[5];
        my $v=join("\t",@f[10..20]);
        $v="$refine_tissue_label\t$v";
        # push @{$hash1{$tissue_ontology_id}},$refine_tissue_label;
        push @{$hash1{$tissue_ontology_id}},$v;
    }
}


foreach my $k(sort keys %hash1){
    my @v1 = @{$hash1{$k}};
    my @uniq_v1 = uniq(@v1);
    # print "$k\t@uniq_v1\n";
    my $num = @uniq_v1;
    if($num >1){
        print "$k\n";
    }
}