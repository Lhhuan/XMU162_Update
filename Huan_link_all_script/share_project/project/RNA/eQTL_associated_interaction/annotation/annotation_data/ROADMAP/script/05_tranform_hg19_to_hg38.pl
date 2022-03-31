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
    print "$marker\n";
    my $fo3 = "${output_dir}/${marker}_sorted_merge.bed.gz";
    my $fo4 = "${output_dir}/hg38/${marker}.bed";
    my $fo5 = "${output_dir}/unmap_hg38/${marker}.bed";   
    my $fo6 = "${output_dir}/hg38/${marker}.bed.gz";
    my $fo7 = "${output_dir}/hg38/${marker}_sorted.bed.gz";
    my $fo8 = "${output_dir}/hg38/${marker}_sorted_merge.bed.gz";
    system "liftOver $fo3 /home/huanhuan/reference/hg19ToHg38.over.chain.gz $fo4 $fo5 ";
    system "gzip $fo4";
    system "zless $fo6 |sort -k1,1 -k2,2n |gzip >$fo7";
    system "bedtools merge -i $fo7 |gzip >$fo8";
}

