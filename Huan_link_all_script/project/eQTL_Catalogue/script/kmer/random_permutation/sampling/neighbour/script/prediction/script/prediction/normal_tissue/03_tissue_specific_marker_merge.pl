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

my $f1 = "./output/02_tissue_level_Marker_source.txt" ;
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^tissue_ontology_id/){
        my @f = split/\t/;
        my $tissue_ontology_id =$f[0];
        my $refine_tissue_label2 =$f[1];
        my $marker =$f[2];
        my $ID =$f[3];
        my $k="$refine_tissue_label2\t$marker";
        push @{$hash1{$k}},$ID;
    }
}

my @histones=("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3");
foreach my $histone(@histones){
    $hash2{$histone}=1
}

foreach my $k(sort keys %hash1){
    my @t=split/\t/,$k;
    my @vs = @{$hash1{$k}};
    my $tissue =$t[0];
    my $marker =$t[1];
    my $outdir= "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/${tissue}/marker"; 
    if (-e $outdir){
        print "$outdir exist";
    }
    else{
        mkpath($outdir);
    }

    #=====
    unless($marker=~/HISTONE_MARK_AND_VARIANT/){
        my $fo2 = "${outdir}/${marker}_merge_pos_info_narrow_peak.bed.gz";
        open my $O2, "| gzip >$fo2" or die $!;
        my $fo3 = "${outdir}/${marker}_merge_pos_info_narrow_peak_signalValue.bed.gz";
        open my $O3, "| gzip >$fo3" or die $!;    
        if(exists $hash2{$marker}){
            foreach my $roadmap_ID(@vs){
                my $f2 = "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}/${roadmap_ID}-${marker}.narrowPeak.gz";         
                open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                while(<$I2>)
                {
                    chomp;
                    my @f =split/\t/;
                    my $chr = $f[0];
                    my $start =$f[1];
                    my $end = $f[2];
                    my $signalValue =$f[6];
                    print $O2 "$chr\t$start\t$end\n";
                    print $O3 "$chr\t$start\t$end\t$signalValue\n";
                }
            }
        }
        else{
            my $dir = $marker;
            $dir =~s/Human_FACTOR/\.\/Human_FACTOR\/human_factor/g;
            $dir =~s/HISTONE_MARK_AND_VARIANT/\.\/HISTONE_MARK_AND_VARIANT\/human_hm/g;
            $dir =~s/Human_CHROMATIN_Accessibility/\.\/Human_CHROMATIN_Accessibility\/human_ca/g;
            foreach my $DCid(@vs){
                my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/$dir/${DCid}_sort_peaks.narrowPeak.bed.gz";         
                open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                while(<$I2>)
                {
                    chomp;
                    my @f =split/\t/;
                    my $chr = $f[0];
                    my $start =$f[1];
                    my $end = $f[2];
                    my $signalValue =$f[6];
                    print $O2 "$chr\t$start\t$end\n";
                    print $O3 "$chr\t$start\t$end\t$signalValue\n";
                }
            }
        }
        close($O2);
        close($O3);
        if(exists $hash2{$marker}){
            my $outdir2 ="$outdir/histone_hg38";
            my $outdir3 ="$outdir/histone_unmap_hg38";
            if (-e $outdir2){
                print "$outdir2 exist";
            }
            else{
                mkdir $outdir2;
            }
            #=====================
            if (-e $outdir3){
                print "$outdir3 exist";
            }
            else{
                mkdir $outdir3;
            }
            #================
            my $fo4="${outdir2}/${marker}_merge_pos_info_narrow_peak.bed";
            my $fo5="${outdir2}/${marker}_merge_pos_info_narrow_peak_signalValue.bed";
            my $fo6="${outdir2}/${marker}_merge_pos_info_narrow_peak_sorted.bed.gz";
            my $fo7="${outdir2}/${marker}_merge_pos_info_narrow_peak_signalValue_sorted.bed.gz";
            system "liftOver $fo2 /home/huanhuan/reference/hg19ToHg38.over.chain.gz $fo4 ${outdir3}/${marker}_merge_pos_info_narrow_peak.bed";
            system "liftOver $fo3 /home/huanhuan/reference/hg19ToHg38.over.chain.gz $fo5 ${outdir3}/${marker}_merge_pos_info_narrow_peak_signalValue.bed";
            system "gzip $fo4";
            system "gzip $fo5";
            system "zless ${fo4}.gz |sort -k1,1 -k2,2n |gzip >$fo6 ";
            system "zless ${fo5}.gz |sort -k1,1 -k2,2n |gzip >$fo7 ";
            system "bedtools merge -i $fo6  |gzip > ${outdir2}/${marker}_merge_pos_info_narrow_peak_sorted_merge.bed.gz"
        }
        else{
            my $fo4="${outdir}/${marker}_merge_pos_info_narrow_peak_sorted.bed.gz";
            my $fo5="${outdir}/${marker}_merge_pos_info_narrow_peak_signalValue_sorted.bed.gz";
            my $fo6="${outdir}/${marker}_merge_pos_info_narrow_peak_sorted_merge.bed.gz";
            system "zless $fo2 |sort -k1,1 -k2,2n |gzip >$fo4";
            system "zless $fo3 |sort -k1,1 -k2,2n |gzip >$fo5";
            system "bedtools merge -i $fo4 |gzip >$fo6";
        }
        print "$k\n";
    }
}
