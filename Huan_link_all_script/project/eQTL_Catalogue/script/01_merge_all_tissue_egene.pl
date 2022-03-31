#利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv" 部分文件，提取pos，p,gene得"../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# my $fo1 = "../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "SNP_chr\tSNP_pos\tPvalue\tegene\n";

my $dir = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data";
my (%hash1,%hash2,%hash3,%hash4);




my $f3 = "/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

while(<$I3>)
{
    chomp;
    unless(/^chr/){
        my @f = split/\t/;
        my $id =$f[0];
        my $maf =$f[1]; 
        my @t=split/:/,$id;
        my $variant =join("\t",@t);
        my $position = join("\t",@t[0..1]);
        $hash3{$variant}=1;
        $hash4{$position}=1;
    }
}



while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $tissue_label = $f[4];
        my $condition_label =$f[5];
        my $quant_method =$f[6];
        my $ftp_path =$f[-2];
        my @t= split/\//,$ftp_path;
        my $filename=$t[-1];
        if($study  =~/\bAlasoo_2018|BLUEPRINT|FUSION|GENCORD|GEUVADIS|HipSci|ROSMAP|Schwartzentruber_2018|TwinsUK|van_de_Bunt_2015|GTEx|Braineac2|iPSCORE|Peng_2018|Steinberg_2020\b/){
            # if($condition_label =~/naive/){
            print "$study\t$filename\n";
            # print "$study\n";
            my $f2 = "${dir}/${study}/${filename}";
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            while(<$I2>)
            {
                chomp;
                unless(/^molecular_trait_id/){
                    my @ss = split/\t/;
                    my $chr = $ss[1];
                    my $pos= $ss[2];
                    my $ref= $ss[3];
                    my $alt= $ss[4];
                    # my $variant =$ss[5];
                    my $pvalue =$ss[8];
                    my $chrpos ="$chr\t$pos";
                    my $variant = "$chrpos\t$ref\t$alt";
                    my $gene_id = $ss[-3];
                    # print "$variant\n$chrpos\n";
                    if(exists $hash3{$variant}){  #test in eur 
                        print $O2 "$chrpos\t$pvalue\t$gene_id\n";

                    }
                    else{
                        if(exists $hash4{$chrpos}){
                            print $O2 "$chrpos\t$pvalue\t$gene_id\n";                              
                        }
                    }
                    
                }
            }
            # }
            # print "$study\n";
            # print $O1 "$_\n";
        }
    }
}

