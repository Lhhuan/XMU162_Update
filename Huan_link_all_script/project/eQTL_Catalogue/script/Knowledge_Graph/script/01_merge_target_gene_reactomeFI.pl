#"../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" 与reactomeFI "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_start_end_all.txt.gz" merge在一起得 ../output/01_hotspot_target_gene_reactomeFI.bed.gz 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f4 = "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_start_end_all.txt.gz"; #
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
my $f5 = "../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"; #
open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件
my $fo1 = "../output/01_hotspot_target_gene_reactomeFI.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "hotspot_chr\thotspot_start\thotspot_end\tegene\tReactomeFI_start_gene\tReactomeFI_end_gene\tReactomeFI_gene\n";
my (%hash1,%hash2);
while(<$I4>)
{
    chomp;
    my @f= split/\t/;
    unless(/^end|start/){
       my $start =$f[0];
       my $end=$f[1];
       push @{$hash1{$start}},$end;
       push @{$hash2{$end}},$start;
    }
}

while(<$I5>)
{
    chomp;
    my @f= split/\t/;
    unless(/^end|start/){
        my $h_chr =$f[0];
        my $h_start= $f[1];
        my $h_end=$f[2];
        my $egene=$f[3];
        # my @all_interaction_gene=();
        my @unique_start=();
        my @unique_end=();
        my (%ha1,%ha2,%ha3);
        if(exists $hash1{$egene}){
            my @ends = @{$hash1{$egene}};
            my @u_ends=grep{++$ha1{$_}<2}@ends;
            push @unique_end,@u_ends;
        }
        if(exists $hash2{$egene}){
            my @starts = @{$hash2{$egene}};
            my @u_starts=grep{++$ha2{$_}<2}@starts;
            push @unique_start,@u_starts;
        }
        my @all_interaction_gene=(@unique_start,@unique_end);
        my @u_all = grep{++$ha3{$_}<2}@all_interaction_gene;
        my $end =join(";",@unique_end);
        my $start =join(";",@unique_start);
        my $all_inter = join(";",@u_all);
        print $O1 "$_\t$start\t$end\t$all_inter\n";
    }
}

