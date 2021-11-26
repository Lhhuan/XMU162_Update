#将每个组织进行显著eQTL-eGene的提取，得"$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz"，排序后和排序后的"../../output/Tissue_total/11_1_extract_max_tissue_share_hotspot_sorted.txt.gz"进行 bedtools intersect,得
#"../../output/Tissue_total/gene/49_share_hotspot_${tissue}_gene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @tissues =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");

my $out_dir ="/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/cis_eQTL_sig/";
my $all_tissue_h = "../../output/cancer_total/share/total/35_1_extract_max_cancer_share_hotspot_sorted.txt.gz";
system "zless ../../output/cancer_total/share/total/35_1_extract_max_cancer_share_hotspot.txt.gz |sort -k1,1 -k2,2n |gzip > $all_tissue_h";
my $anno_out_dir = "../../output/cancer_total/share/total/gene";
my $pm = Parallel::ForkManager->new(10);
foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    print "$tissue\tstart\n";
    my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/cis_eQTLs_all_re.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $fo1 = "$out_dir/${tissue}_cis_sig_eQTL_egene.txt.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    # print $O1 "SNP_chr\tSNP_pos\tPvalue\n";


    my %hash1;
    while(<$I1>)
    {
        chomp;
        unless(/^cancer_type/){
            my @f = split/\t/;
            my $cancer_type=$f[0];
            my $SNP_chr =$f[2];
            my $SNP_pos =$f[3];
            my $gene_id1 = $f[5];
            my $Pvalue =$f[-1];
            my $start = $SNP_pos;
            my $end = $SNP_pos +1;
            my $gene_id = $gene_id1;
            # $gene_id =~ s/\..*+//g;
            # print "$gene_id1\t$gene_id\n";
            if($Pvalue < 5E-8 && $cancer_type eq $tissue  ){
                my $output = "$SNP_chr\t$start\t$end\t$gene_id";
                unless(exists $hash1{$output}){
                    $hash1{$output}=1;
                    print $O1 "$output\n";
                }
            }
        }
    }
    close($O1);
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip > $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz";
    system "bedtools intersect -a $all_tissue_h -b $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz -wo |cut -f1-3,7|sort -u |gzip > $anno_out_dir/19_share_hotspot_${tissue}_gene.txt.gz";
    # system "bedtools intersect -a $all_tissue_h -b $out_dir/${tissue}_cis_sig_eQTL_egene_sorted.txt.gz -wo |sort -u |gzip > $anno_out_dir/49_share_hotspot_${tissue}_gene.txt.gz" ;
    print "$tissue\tend\n";
    $pm->finish;  #多线程结束
}



