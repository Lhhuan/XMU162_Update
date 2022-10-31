
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

my $f1 = "../output/01_hotspot_idx.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/01_gene_idx.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $f3 = "../../../../../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz";
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
my $f4 = "../../../../../output/nodes_annotation/hotspot_annotation.txt.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

my $fo1 = "../output/02_hotspot_egene_idx.txt.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Hotspot_idx\tEgene_idx\n";
my $fo2 = "../output/02_egene_interaction_idx.txt.gz"; #
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "Egene_idx\tReactomeFI_idx\n";
my $fo3 = "../output/02_egene_co-expression_idx.txt.gz"; #
open my $O3, "| gzip >$fo3" or die $!;
print $O3 "Egene_idx\tco_expression_idx\n";
my $fo4 = "../output/02_hotspot_egene_reactomeFI_co-expression_idx.bed.gz"; #
open my $O4, "| gzip >$fo4" or die $!;
print $O4 "Hotspot_idx\tegene_idx\tReactomeFI_gene_idx\tegene_co_expression_gene_idx\n";

my $fo5 = "../output/02_hotspot_anno_idx.txt.gz"; #
open my $O5, "| gzip >$fo5" or die $!;
print $O5 "Hotspot_idx\tH3K27ac\tH3K4me1\tH3K4me3\tH3K9ac\tH3K36me3\tH3K27me3\tH3K9me3\tCTCF\tCHROMATIN_Accessibility\tTFBS\tEnhancer\n";

my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    unless(/^Hotspot/){
        my @f =split/\t/;
        my $hotspot  =$f[0];
        my $idx=$f[1];
        $hash1{$hotspot}=$idx;
    }
}

while(<$I2>)
{
    chomp;
    unless(/^Gene/){
        my @f =split/\t/;
        my $Gene  =$f[0];
        my $idx=$f[1];
        $hash1{$Gene}=$idx;
    }
}

while(<$I3>)
{
    chomp;
    unless(/^hotspot_chr/){
        my @f =split/\t/;
        my $hotspot=join("_",@f[0..2]);
        if(exists $hash1{$hotspot}){
            my $hotspot_idx = $hash1{$hotspot};
            my $egene =$f[3];
            my $egene_idx =$hash1{$egene};
            my $ReactomeFI_gene =$f[6];
            my $egene_co_expression_gene =$f[7];
            my @interaction=();
            my @co=();
            print $O1 "$hotspot_idx\t$egene_idx\n";
            unless($ReactomeFI_gene =~/NA/){
                my @ts= split/;/,$ReactomeFI_gene;
                foreach my $t(@ts){
                    my $gene_idx=$hash1{$t};
                    push @interaction,$gene_idx;
                }
                my $ReactomeFI_gene_idx =join(";",@interaction);
                my $output2= "$egene_idx\t$ReactomeFI_gene_idx";
                unless(exists $hash2{$output2}){
                    $hash2{$output2}=1;
                    print $O2 "$output2\n";
                }
            }
            unless($egene_co_expression_gene =~/NA/){
                my @ts= split/;/,$egene_co_expression_gene;
                foreach my $t(@ts){
                    my $gene_idx=$hash1{$t};
                    push @co,$gene_idx;
                }            
                my $egene_co_expression_gene_idx =join(";",@co);
                my $output3 ="$egene_idx\t$egene_co_expression_gene_idx";
                unless(exists $hash3{$output3}){
                    $hash3{$output3}=1;
                    print $O3 "$output3\n";
                }
                # print $O3 "$egene_idx\t$egene_co_expression_gene_idx\n";
            }
            my $ReactomeFI_gene_idx =join(";",@interaction);
            my $egene_co_expression_gene_idx =join(";",@co);
            print $O4 "$hotspot_idx\t$egene_idx\t$ReactomeFI_gene_idx\t$egene_co_expression_gene_idx\n";
        }
    }
}

while(<$I4>)
{
    chomp;
    unless(/^Hotspot/){
        my @f=split/\t/;
        my $hotspot=$f[0];
        my $marker=join("\t",@f[1..$#f]);
        if(exists $hash1{$hotspot}){        
            my $hotspot_idx =$hash1{$hotspot};
            print $O5 "$hotspot_idx\t$marker\n";
        }
    }
}