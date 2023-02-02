#根据"/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv" merge数据

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::MoreUtils ':all';

my $f1 = "/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/02_hotspots_tissue_type_annotation_cistromeDB_info.tsv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 



my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $marker = $f[2];
        my $DCid =$f[3];
        my $cell_line =$f[4];
        # my $v="$DCid\t$cell_line";
        push @{$hash1{$marker}},$DCid;
    }
}

my %hash5;
for (my $i=1;$i<23;$i++){
    my $k = "chr${i}";
    $hash5{$k}=1;
}

my $fo4 = "./normal_cell/052_normal_hotspot_narrowpeak_id.txt";
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";

print $O4 "marker\tDCid\n";


foreach my $k(sort keys %hash1){
    print "$k\tstart\n";
    my $dir = $k;
    $dir =~s/Human_FACTOR/\.\/Human_FACTOR\/human_factor/g;
    $dir =~s/HISTONE_MARK_AND_VARIANT/\.\/HISTONE_MARK_AND_VARIANT\/human_hm/g;
    $dir =~s/Human_CHROMATIN_Accessibility/\.\/Human_CHROMATIN_Accessibility\/human_ca/g;
    my @vs =@{$hash1{$k}};
    my @uniq_vs=uniq(@vs);
    my $output_dir = "./normal_cell/${k}";
    my $fo1 = "${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "${output_dir}/merge_pos_info_narrow_peak.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    my $fo3 = "${output_dir}/merge_pos_info_narrow_peak_signalValue.bed.gz";
    open my $O3, "| gzip >$fo3" or die $!;
    my(%hash2,%hash3,%hash4);
    foreach my $v(@uniq_vs){
        $hash2{$v}=1;
    }
    opendir (DIR, $dir) or die "can't open the directory!";
    my @files = readdir DIR;
    foreach my $file(@files){
        if ( $file =~ /[a-z]/) {
            my @t = split/\_/,$file;
            my $id = $t[0];
            # print "$id\n";
            if (exists $hash2{$id}){
                if ($file =~/narrowPeak/){
                    print $O4 "$k\t$id\n";
                    my $f2 = "$dir/$file";
                    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                    while(<$I2>)
                    {
                        chomp;
                        my @f = split/\s+/;
                        my $chr = $f[0];
                        my $start =$f[1];
                        my $end =$f[2];
                        my $name = $f[3];
                        my $signalValue =$f[6];
                        if (exists $hash5{$chr}){
                            my $output1 = "$chr\t$start\t$end\t$id";
                            my $output2 = "$chr\t$start\t$end";
                            my $output3 = "$chr\t$start\t$end\t$signalValue";
                            print $O3 "$output3\n";
                            unless(exists $hash3{$output1}){
                                $hash3{$output1}=1;
                                print $O1 "$output1\n";
                            }
                            unless(exists $hash4{$output2}){
                                $hash4{$output2}=1;
                                print $O2 "$output2\n";
                            }
                        }
                    }
                }
            }
        }
    }
    close($O1);
    close($O3);
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ${output_dir}/merge_pos_info_sample_narrow_peak_sorted.bed.gz";
    system "zless $fo3 |sort -k1,1 -k2,2n |gzip > ${output_dir}/merge_pos_info_narrow_peak_signalValue_sorted.bed.gz";
    print "$k\tfinish\n";
}
