#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager;   

my (%hash1,%hash2);
for(my $i=1; $i<23;$i++){
    $hash1{$i}=1;
}
    
my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Lung_cis_eQTL_1kg_Completion.txt.gz";
# my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/123.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my $fo1 = "Lung_cis_eQTL_1kg_Completion_chr_1_22_unique_pos.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "SNP_chr\tSNP_pos\n";
 
while(<$I2>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $SNP_chr =$f[0];
        my $SNP_pos =$f[1];
        my $Pvalue =$f[2];
        if (exists $hash1{$SNP_chr}){
            my $output = "$SNP_chr\t$SNP_pos";
            unless (exists $hash2{$output}){
                $hash2{$output}=1;
                print $O1 "$output\n";
            }
        }
        
    }
} 