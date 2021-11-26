#将不同sample的相同mark进行合并，得 "/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_merge.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;


my @markers = ("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
# my @markers = ("H3K4me3");
# my @markers = ("H3K4me3","H3K9ac");

my $output_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/merge";
foreach my $marker(@markers){
    my $fo1 = "${output_dir}/${marker}.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "${output_dir}/${marker}_sorted.bed.gz";
    my $fo3 = "${output_dir}/${marker}_sorted_merge.bed.gz";   


    my $f1 = "unique_roadmap_id.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    while(<$I1>)
    {
        chomp;
        my $roadmap_ID = $_;
        my $sample_marker_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}";
        my $f2 = "$sample_marker_dir/${roadmap_ID}-${marker}.narrowPeak.gz";
        if(-e $f2){
            print "$roadmap_ID\t$marker\n";
            # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            while(<$I2>)
            {
                chomp;
                my @f =split/\t/;
                my $chr = $f[0];
                my $start =$f[1];
                my $end = $f[2];
                print $O1 "$chr\t$start\t$end\n";
            }
        }
    }
    close($f1);
    close($O1);
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip >$fo2";
    system "bedtools merge -i $fo2 |gzip >$fo3";
}

