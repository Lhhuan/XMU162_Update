#将"/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_sample.bed.gz" 从 hg19转换到hg38,得"/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/${marker}_sorted_sample.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
# use Parallel::ForkManager;


my @markers = ("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
# my @markers = ("H3K4me3");
# my @markers = ("H3K4me3","H3K9ac");
my $sample = "E062";
my $tissue = "Blood";
my $input_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/${sample}";
my $output_dir ="/share/data0/GTEx/annotation/ROADMAP/sample/tissue_specific/${tissue}";
foreach my $marker(@markers){
    print "$marker\n";
    my $input = "${input_dir}/${sample}-${marker}.narrowPeak.gz";
    my $fo3 = "${output_dir}/${marker}_signalValue_sorted.bed.gz";
    my $fo4 = "${output_dir}/hg38/${marker}_signalValue.bed";
    my $fo5 = "${output_dir}/unmap_hg38/${marker}_signalValue.bed";   
    my $fo6 = "${output_dir}/hg38/${marker}_signalValue.bed.gz";
    my $fo7 = "${output_dir}/hg38/${marker}_signalValue_sorted.bed.gz";
    my $fo8 = "${output_dir}/hg38/${marker}_sorted.bed.gz";
    my $fo9 = "${output_dir}/hg38/${marker}_sorted_merge.bed.gz";
    system "zless $input |cut -f1-3,7 |sort -k1,1 -k2,2n|gzip >$fo3";
    system "liftOver $fo3 /home/huanhuan/reference/hg19ToHg38.over.chain.gz $fo4 $fo5 ";
    system "gzip $fo4";
    system "zless $fo6 |sort -k1,1 -k2,2n |gzip >$fo7";
    system "zless $fo7 |cut -f1-3 | gzip >$fo8";
    system "bedtools merge -i $fo8 |gzip >$fo9"
}

