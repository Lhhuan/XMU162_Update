#根据"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" 的 tissue_ontology_id在"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_roadmap_info.tsv" 和"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv"提取tissue 的annotation info,得./output/02_tissue_level_Marker_source.txt, refine marker source 得"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info_tissue_refine.tsv"
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
my $f2 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_roadmap_info.tsv" ;
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $f3 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv" ;
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
my $fo1 = "./output/02_tissue_level_Marker_source.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info_tissue_refine.tsv";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
print $O1 "tissue_ontology_id\trefine_tissue_label2\tmarker\tID\n";


my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $before = join("\t",@f[0..5]);
    my $later =join("\t",@f[6..$#f]);
    if(/^study/){
        print $O2 "$before\trefine_tissue_label2\t$later\n";
    }
    else{
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $tissue_ontology_id =$f[2];
        my $refine_tissue_label =$f[5];
        my $refine_tissue_label2 =$refine_tissue_label; 
        $refine_tissue_label2 =~s/\(//g;
        $refine_tissue_label2 =~s/\)//g;
        my $k="$tissue_ontology_id\t$refine_tissue_label2";
        my $v= "$study\t$qtl_group";
        push @{$hash1{$k}},$v;

        print $O2 "$before\t$refine_tissue_label2\t$later\n";
    }
}

while(<$I2>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $roadmap_id =$f[2];
        my $marker =$f[3];
        my $k="$study\t$qtl_group";
        my $v= "$marker\t$roadmap_id";
        push @{$hash2{$k}},$v;
    }
}

while(<$I3>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $marker =$f[2];
        my $DCid =$f[3];
        my $k="$study\t$qtl_group";
        my $v= "$marker\t$DCid";
        push @{$hash2{$k}},$v;
    }
}

my %hash3;
foreach my $k1(sort keys %hash1){
    my @v1s = @{$hash1{$k1}};
    foreach my $v1(@v1s){
        if (exists $hash2{$v1}){
            my @v2s = @{$hash2{$v1}};
            my @uniq_v2s = uniq(@v2s);
            foreach my $v2(@uniq_v2s){
                my $output = "$k1\t$v2";
                unless(exists $hash3{$output}){
                    $hash3{$output}=1;
                    print $O1 "$output\n";
                }
            }
        }
    }
}