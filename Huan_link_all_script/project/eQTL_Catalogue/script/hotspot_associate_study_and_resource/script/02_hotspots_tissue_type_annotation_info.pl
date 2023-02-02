#为"../output/01_adjust_tissue_need_study_for_hotspot_annotation_cistromeDB.tsv"的roadmap和cistromedb提取信息分别得"../output/02_hotspots_tissue_type_annotation_roadmap_info.tsv"，"../output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv",并得去掉引号的文件"../output/02_adjust_tissue_need_study_for_hotspot_annotation_cistromeDB_refine.tsv"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/051_normal_chrom_marker_info_refine.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 = "../output/01_adjust_tissue_need_study_for_hotspot_annotation_cistromeDB.tsv" ;
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $fo1 = "../output/02_hotspots_tissue_type_annotation_roadmap_info.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "../output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 = "../output/02_adjust_tissue_need_study_for_hotspot_annotation_cistromeDB_refine.tsv";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";

print $O1 "study\tqtl_group\troadmap_id\tMarker\n";
print $O2 "study\tqtl_group\tMarker\tDCid\tcell_line\tRefine_tissue_type (Refine_Cell_type)\n";


my (%hash1,%hash2,%hash3);
while(<$I1>)
{
    chomp;
    unless(/^Marker/){
        my @f=split/\t/;
        my $Marker =$f[0];
        my $DCid =$f[1];
        my $cell_line=$f[-3];
        my $Refine_Cell_type =$f[-2];
        my $Refine_tissue_type =$f[-1];
        my $v= "$Marker\t$DCid\t$cell_line";
        unless($Refine_Cell_type =~/NA/){
            my $k ="$Refine_tissue_type ($Refine_Cell_type)";
            $k =~ s/"//g;
            push @{$hash1{$k}},$v;
            # print "$k\n";
        }
    }
}


while(<$I2>)
{
    chomp;
    if(/^study/){
        print $O3 "$_\n";
    }
    else{
        my $line = $_;
        $line =~s/"//g;
        print $O3 "$line\n";
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $roadmap= $f[-2];
        $roadmap =~ s/"//g;
        my $cistromeDB =$f[-1];
        $cistromeDB =~ s/"//g;
        my $k = "$study\t$qtl_group";
        # my $v1 = join("\t",@f[0..$#f]);
        unless($roadmap =~/NA/){
            my @roadmaps = split/\,/,$roadmap;
            foreach my $histone(@roadmaps){
                # print "$k\t$histone\n"
                my @t =split/\(/,$histone;
                my $id=$t[1];
                $id =~s/\)//g;
                # print "$id\n";
                push @{$hash2{$k}},$id; 
            }
        }
        unless($cistromeDB =~/NA/){
            my @cistromeDBs = split/\,/,$cistromeDB;
            foreach my $tfbs(@cistromeDBs){
                if(exists $hash1{$tfbs}){
                    my @infos = @{$hash1{$tfbs}};
                    foreach my $info(@infos){
                        my $v3 = "$info\t$tfbs";
                        push @{$hash3{$k}},$v3;
                        # print "$tfbs\t$info\n";
                    }
                }
            }
        }
    }
}
my @markers = ("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");

foreach my $k(sort keys %hash2){
    my @vs = @{$hash2{$k}};
    foreach my $roadmap_ID(@vs){
        foreach my $marker(@markers){
            my $marker_file = "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}/${roadmap_ID}-${marker}.narrowPeak.gz";
            if(-e $marker_file){
                print $O1 "$k\t$roadmap_ID\t$marker\n";
            }
            
            # print "$roadmap_ID\t$marker\n";
        }
    }
}

foreach my $k(sort keys %hash3){
    my @vs = @{$hash3{$k}};
    foreach my $v(@vs){
        print $O2 "$k\t$v\n";
    }   
}