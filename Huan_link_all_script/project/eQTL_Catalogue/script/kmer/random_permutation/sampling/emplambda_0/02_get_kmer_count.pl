#将$input_file 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $out_dir = "/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/emp0/";
my $pm = Parallel::ForkManager->new(30);
for (my $i=1;$i<300;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $bed_input = "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/sorted_${i}_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
    my $output_fa = "${out_dir}/${i}_emp0_hotspot_random.fa";
    my $command1 = " bedtools getfasta -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed  $bed_input -fo $output_fa";
    system $command1;

    chdir $out_dir;
    $ENV{'i'} = $i;
    $ENV{'output_fa'} = $output_fa;
    my $command2 = "bash /home/huanhuan/project/eQTL_Catalogue/script/kmer/get_kmer.sh";
    system $command2;
    print "$i\n";
    $pm->finish;  #多线程结束
}

