#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
# use Parallel::ForkManager;

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
for(my $i=1;$i<8;$i++){
    my $f1 = "../output/figure/09_chr/11_chr1_louvain_pca5_k500.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    my $out_dir1 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200";
    if(-e $out_dir1){
        print "${out_dir1}\texist\n";
    }
    else{
        system "mkdir -p $out_dir1";
    }
    my $fo1 = "${out_dir1}/cluster_${i}.bed";
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
    my $output_dir = "${out_dir1}/${i}";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    system "findMotifsGenome.pl $fo1 hg38 $output_dir -size 200";
    print "$i\tfinish\n";
}