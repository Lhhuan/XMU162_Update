#"../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 与reactomeFI "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part1_start_end.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part2_end_start.txt.gz","/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_start_end_all.txt.gz" merge在一起得 ../output/01_hotspot_target_gene_reactomeFI.bed.gz 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

# my $f1 = "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part1_start_end.txt.gz";
# # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# my $f2 = "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part2_end_start.txt.gz"; #
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
# my $f3 = "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz"; #
# open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
my $f4 = "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_start_end_all.txt.gz"; #
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
# my $f5 = "../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"; #
my $f5 = "1234.txt.gz"; 
open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件
my $fo1 = "../output/01_ori_hotspot_target_gene_reactomeFI.bed.gz"; #
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
       if(exists $hash1{$egene}){
           my @ends = @{$hash1{$egene}};
           my @starts = @{$hash2{$egene}};
           my @all_interaction_gene=(@starts,@ends);
           my (%ha1,%ha2,%ha3);
           my @u_ends=grep{++$ha1{$_}<2}@ends;
           my @u_starts=grep{++$ha2{$_}<2}@starts;
           my @u_all = grep{++$ha3{$_}<2}@all_interaction_gene;
           my $end =join(";",@u_ends);
           my $start =join(";",@u_starts);
           my $all_inter = join(";",@u_all);
           print $O1 "$_\t$start\t$end\t$all_inter\n";
       }
        if(exists $hash1{$egene}){
            print "$egene\n";
        }
        else{
            print "$egene\tnot_exists\n";
        }
    }
}

