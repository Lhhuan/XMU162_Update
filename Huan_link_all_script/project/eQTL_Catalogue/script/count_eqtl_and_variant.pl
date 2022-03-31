#合并部分数据集位点并利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz"筛选落在eur maf >0.05的eqtl,同时利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.01.txt.gz" 补全得"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz"
#并得进入筛选的数据集信息"../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "../output_morethan0.01/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output_morethan0.01/01_merge_all_tissue_cis_eQTL_1kg_Completion_unique_eQTL_pos.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output_morethan0.01/01_merge_all_tissue_cis_eQTL_1kg_Completion_unique_variant_pos.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
# print $O2 "SNP_chr\tSNP_pos\tPvalue\n";

my $dir = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data";
my (%hash1,%hash2,%hash3,%hash4);





while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f=split/\t/;
        my $chr = $f[0];
        my $SNP_pos =$f[1];
        my $Pvalue =$f[2];
        my $output= "$chr\t$SNP_pos";
        if($Pvalue <=5e-8){
            unless(exists $hash1{$output}){
                $hash1{$output}=1;
                print $O1 "$output\n";
            }
        }
        unless(exists $hash2{$output}){
            $hash2{$output}=1;
            print $O2 "$output\n";
        }        
    }
}

