#将$input_file 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $out_dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/";
my $pm = Parallel::ForkManager->new(20);
for (my $i=301;$i<401;$i++){
    # my $i=100;
    my $pid = $pm->start and next; #开始多线程
    my $bed_input = "${out_dir}/${i}_10_fold_neighbour_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
    my $sorted_bed = "${out_dir}/${i}_sorted_10_fold_neighbour_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
    my $command01 = "zless ${bed_input} |sort -k1,1 -k2,2n |gzip >${sorted_bed} ";
    system $command01;
    my $output_fa = "${out_dir}/${i}_emp10_fold_hotspot_random.fa";
    my $command1 = " bedtools getfasta -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed  $sorted_bed -fo $output_fa";
    system $command1;

    chdir $out_dir;
    $ENV{'i'} = $i;
    $ENV{'output_fa'} = $output_fa;
    my $command2 = "bash /home/huanhuan/project/eQTL_Catalogue/script/kmer/get_kmer.sh";
    system $command2;
    print "$i\n";
    $pm->finish;  #多线程结束
}

