#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 

my $out_dir1 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/GC_content/";
if(-e $out_dir1){
    print "${out_dir1}\texist\n";
}
else{
    system "mkdir -p $out_dir1";
}


my $fo2 = "${out_dir1}/gc_content.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "hotspot_chr\thotspot_start\thotspot_end\tgc_content\tcluster\n";



for(my $i=1;$i<7;$i++){
    my $f1 = "../../../output/figure/whole_genome/11_whole_genome_leiden_pca3_k50_resolution1e-04.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    my $fo1 = "${out_dir1}/${i}_cluster.bed";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    while(<$I1>)
    {
        chomp;
        unless(/^PC_1/){
            my @f= split/\t/;
            my $cluster=$f[-2];
            my $hotspot=$f[-1];
            my @t=split/\:/,$hotspot;
            my $chr = $t[0];
            my @ss=split/-/,$t[1];
            my $start = $ss[0];
            my $end =$ss[1];
            # print "$cluster\n";
            if($cluster =~/$i/){
                my $output = "$chr\t$start\t$end";
                print $O1 "$output\n";
            }
        }
    } 
    close($O1);
    my $output_dir = ${out_dir1};
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    my $gc_file = "${out_dir1}/${i}_cluster_gc_content.bed.gz";
    system "bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed $fo1 |gzip >$gc_file";
    print "$i\tfinish\n";
    my $f2=$gc_file;
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); 
    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f= split/\t/;
            my $hotspot=join("\t",@f[0..2]);
            my $gc_count=$f[4];
            print $O2 "$hotspot\t$gc_count\t$i\n";
        }
    } 
}