#利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv" 部分文件，过滤p<0.05得"../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05.txt.gz"
#过滤p<5e-8得"../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "SNP_chr\tSNP_start\tSNP_start\tegene\tPvalue\n";

my (%hash1,%hash2,%hash3,%hash4);

while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $SNP_chr =$f[0];
        my $SNP_pos =$f[1];
        my $Pvalue =$f[2];
        my $egene =$f[3];
        my $tissue=$f[4];
        my $chr ="chr${SNP_chr}";
        my $start = $SNP_pos-1+1;
        my $end =$SNP_pos+1;
        my $output = "$chr\t$start\t$end\t$egene\t$Pvalue\t$tissue";
        if($Pvalue <0.05){
            unless(exists $hash1{$output}){
                $hash1{$output}=1;
                print $O1 "$output\n";
            }
            if($Pvalue <5e-8){
                # print "$_\n";
                unless(exists $hash2{$output}){
                    $hash2{$output}=1;
                    print $O2 "$output\n";
                }      
            }
        }

    }
}


