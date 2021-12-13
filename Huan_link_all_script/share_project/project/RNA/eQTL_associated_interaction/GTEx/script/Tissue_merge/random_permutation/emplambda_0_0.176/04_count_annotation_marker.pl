#对/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/$type/$cutoff2/*的marker进行count,得"$out_dir/${type}_${cutoff2}_histone_marker.txt.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Tissue_merge";
my $type = "original_random";
# my $type = "background_original_random";
my $cutoff2 = "0_0.176";

my $input_dir = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/interval_18_random/original_random/extend/emp${cutoff}/${cutoff2}";

# my @markers = ("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3");
my @markers = ("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","CHROMATIN_Accessibility","TFBS","CTCF");
my $out_dir=  "../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/enrichment/";



my $fo1 = "$out_dir/${type}_${cutoff2}_histone_marker.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Marker\tcount\trandom_number\n";




# my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
for (my $i=1;$i<10001;$i++){
# for (my $i=500;$i<1001;$i++){
    # my $pid = $pm->start and next; #开始多线程
    foreach my $mark (@markers){
        my $input_file_base_name = "${mark}_${i}_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz";
        my $input_file = "$input_dir/$input_file_base_name";
        #----------------------------
        my @arg1s = stat ($input_file);
        my $mark_size = $arg1s[7]; 
        my $command_mark = "zless $input_file |cut -f1,2,3| sort -u | wc -l" ;
        my $command_mark_line_count = wc($command_mark);
        # print "$mark\t$command_mark_line_count\t$i\n";
        print $O1 "$mark\t$command_mark_line_count\t$i\n";
    }
    print "$i\n";
}
# print "$type\tend\n";




sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}