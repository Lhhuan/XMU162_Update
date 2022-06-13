#../data/endothelial_cell_of_hepatic_sinusoid/ENCFF879NDS.bed.gz 转成hg38,再合并到一起得"../output/hg38/01_merge_TAD_sample.bed.gz"，排序得../output/hg38/01_merge_TAD_sample_sorted.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use File::Basename;
my (%hash1,%hash2);

my $fo1 = "../output/hg38/01_merge_TAD_sample.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;


my @files= ("../data/endothelial_cell_of_hepatic_sinusoid/ENCFF879NDS.bed.gz","../data/endometrial_microvascular_endothelial_cells/ENCFF633ORE.bed.gz","../data/astrocyte_of_the_cerebellum/ENCFF306YQN.bed.gz","../data/astrocyte_of_the_spinal_cord/ENCFF444XGA.bed.gz");

foreach my $file(@files){
    my @t=split/\//,$file;
    my $tissue= $t[2];
    my $file_name = basename($file);
    my $dir = dirname($file);
    system "mkdir -p $dir/hg38";
    my $f3="./$dir/hg38/01_${tissue}.bed";
    system "liftOver $file  /home/huanhuan/reference/hg19ToHg38.over.chain.gz $f3 ./$dir/01_unmap_${tissue}.bed";       
    open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
    while(<$I3>)
    {
        chomp;
        my @f = split/\t/;
        my $pos= join("\t",@f[0..2]);
        print $O1 "$pos\t$tissue\n";
        
    }   
}
close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ../output/hg38/01_merge_TAD_sample_sorted.bed.gz";