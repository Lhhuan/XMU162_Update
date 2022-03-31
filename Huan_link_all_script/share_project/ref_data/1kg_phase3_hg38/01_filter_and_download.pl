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



my $f1 ="phased-manifest_July2021.tsv";; 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "need_download_file.tsv";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my @all_link=();
while(<$I1>){
    chomp;
    my @f = split/\t/;
    my $file_name= $f[0];
    # if($file_name=~/variants.vcf/){
    my $link = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${file_name}";
    print $O1 "$link\n";
    push @all_link,$link;
    # }
}

my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
foreach my $link(@all_link){
    my $pid = $pm->start and next; #开始多线程
    system "wget -c $link";
    $pm->finish;  #多线程结束
    
}