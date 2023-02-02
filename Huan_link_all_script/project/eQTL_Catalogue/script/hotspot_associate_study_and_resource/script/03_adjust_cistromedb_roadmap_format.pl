#将"../output/02_hotspots_tissue_type_annotation_roadmap_info.tsv"，"../output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv",调整合并得"../output/03_hotspots_tissue_type_annotation_roadmap_cistromeDB_refine.tsv"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::MoreUtils ':all';
my $f1 = "../output/02_hotspots_tissue_type_annotation_roadmap_info.tsv" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 = "../output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv" ;
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $fo1 = "../output/03_hotspots_tissue_type_annotation_roadmap_cistromeDB_refine.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";


print $O1 "study\tqtl_group\tMarker\tCell_or_RoadmapID\n";



my (%hash1,%hash2,%hash3);
while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f=split/\t/;
        my $study =$f[0];
        my $qtl_group =$f[1];
        my $roadmap_id =$f[2];
        my $Marker =$f[3];
        my $k = "$study\t$qtl_group\t$Marker";
        push @{$hash1{$k}},$roadmap_id;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^study/){
        my @f=split/\t/;
        my $study =$f[0];
        my $qtl_group =$f[1];
        my $Marker =$f[2];
        my $cell_line =$f[4];

        my $k = "$study\t$qtl_group\t$Marker";
        push @{$hash1{$k}},$cell_line;
    }
}

foreach my $k(sort keys %hash1){
    my @vs = @{$hash1{$k}};
    my @uniq_vs =uniq(@vs);
    my $v=join(";",@uniq_vs);
    print $O1 "$k\t$v\n";
}


