#将所有tissue "${dir}/${tissue}${suffix}"合并得gene文件"../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz";
#合并位点并用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全得"../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
use Env qw(PATH);



my $f1 ="tabix_ftp_paths.tsv";; 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "need_download_tabix_ftp_paths.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
while(<$I1>){
    chomp;
    if(/^study/){
        print $O1 "$_\trefine_tissue_label\n";
    }
    else{
        my @f = split/\t/;
        my $study= $f[0];
        my $tissue_label = $f[4];
        my $condition_label =$f[5];
        my $quant_method =$f[6];
        my $ftp_path =$f[-1];
        # if($condition_label =~/naive/){
            # print "$tissue_label\n";
            # unless($tissue_label =~/\bcell|monocyte|fibroblast|iPSC|macrophage|microglia|monocyte|neutrophil|platelet|Treg\b/){
                if($quant_method=~/ge/){
                    my $tissue_label1 =$tissue_label;
                    $tissue_label1 =~ s/\s+/_/g;
                    print "$study\t$tissue_label1\n";
                    print $O1 "$_\t$tissue_label1\n";
                    $ENV{'study'}  = $study; #设置环境变量
                    $ENV{'ftp_path'} = $ftp_path;
                    system "bash download.sh";
                }
            # }
        # }
    }
}


