#根据"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/052_normal_hotspot_narrowpeak_id.txt" 过滤 "../output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv",得"../output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/052_normal_hotspot_narrowpeak_id.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $f2 = "../output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv" ;
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $fo1 = "../output/021_hotspots_tissue_type_annotation_cistromeDB_info_filter_narrowpeak.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";



my (%hash1,%hash2,%hash3);
while(<$I1>)
{
    chomp;
    unless(/^marker/){
        my @f=split/\t/;
        my $Marker =$f[0];
        my $DCid =$f[1];
        my $k = "$Marker\t$DCid";
        $hash1{$k}=1;
    }
}


while(<$I2>)
{
    chomp;
    if(/^study/){
        print $O1 "$_\n";
    }
    else{
        my @f=split/\t/;
        my $Marker =$f[2];
        my $DCid =$f[3];
        my $k="$Marker\t$DCid";
        if(exists $hash1{$k}){
            print $O1 "$_\n";
        }
    }
}