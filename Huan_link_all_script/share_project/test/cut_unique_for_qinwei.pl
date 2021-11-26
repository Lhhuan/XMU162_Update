#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/qinwei/66panel_anno.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");


my $fo1 = "intersect_unique.bed";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;

my %hash1;


while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $chr = $f[0];
    my $type = $f[5];
    my $start = $f[1];
    my $end  =$f[2];
    my $info =$f[11];
    if ($type =~/exon/){
        $info =~ s/\s+//g;
        my @ts = split/;/,$info;
        my @ele= ();
        foreach my $t(@ts){
            if ($t =~ /^gene_id/ || $t =~/^gene_name/){
                
                $t =~ s/gene_id|"|gene_name//g;
                # print "$t\n";
                push @ele, $t;
            }
        }
        my $final_ele = join ("\t",@ele);
        my $output ="$chr\t$start\t$end\t$final_ele\t$type";
        unless(exists $hash1{$output}){
            $hash1{$output} =1;
            print $O1 "$output\n";
        }
    }

}
