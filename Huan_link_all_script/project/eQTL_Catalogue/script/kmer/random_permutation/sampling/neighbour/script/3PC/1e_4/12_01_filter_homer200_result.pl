#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 

my $q = 0.05;
my $out_dir1 = "/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/homer_200/";
my $fo1 = "${out_dir1}/homer200_result_qvalue${q}.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "${out_dir1}/homer200_result_cluster_specific_motif_qvalue${q}.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
print $O1 "Motif_name\tConsensus\tpvalue\tqvalue\tCluster\n";
print $O2 "motif_name\tCluster\n";

my %hash1;
foreach my $i(1..6){
    my $f1 = "${out_dir1}/${i}/knownResults.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    while(<$I1>)
    {
        chomp;
        unless(/^Motif/){
            my @f= split/\t/;
            my $Motif_name =$f[0];
            # $Motif_name =~ s/\(.*//g;
            $Motif_name =~ s/\/.*//g;
            # print "$Motif_name\n";
            my $Consensus =$f[1];
            my $pvalue= $f[2];
            my $qvalue= $f[4];
            if($qvalue <$q){
                print $O1 "$Motif_name\t$Consensus\t$pvalue\t$qvalue\t$i\n";
                push @{$hash1{$Motif_name}},$i
            }
        }
    } 
    close($I1);
    print "$i\tfinish\n";
}

foreach my $k(keys %hash1){
    my @vs = @{$hash1{$k}};
    my @uni=uniq(@vs);
    my $num=@uni;
    if($num <2){
        print $O2 "$k\t$uni[0]\n";
    }
}