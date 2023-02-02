#根据"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/04_merge_hotspot_and_annotation_info.tsv" 的 tissue_ontology_id在"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_roadmap_info.tsv" 和"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv"提取tissue 的annotation info,得./output/02_tissue_level_Marker_source.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::MoreUtils ':all';
use File::Path;
# my $dir ="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/aaaa/bbbb";
# mkpath($dir);
my $f1 = "./output/02_tissue_level_Marker_source.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^tissue_ontology_id/){
        my @f = split/\t/;
        my $refine_tissue_label2 =$f[1];
        $hash1{$refine_tissue_label2}=1;
    }
}

foreach my $k(sort keys %hash1){
    print "$k\n";
}

# my $f1 = "1234.tsv" ;
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# # my $fo1 = "./output/01_adjust_tissue_need_study_for_hotspot.tsv";
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

# if(-z $f1){
#     print "f1 is null\n";
# }
# # esle



# # my %hash1;

# while(<$I1>)
# {
#     chomp;
#     print "HHHH\n";
#     # unless(/^study/){
#         my @f = split/\t/;
#         my $chr = $f[0];
#         my $start =$f[1];
#         my $end =$f[2];
#         my $name = $f[3];
#         my $signalValue =$f[6];
#     # }
# }
