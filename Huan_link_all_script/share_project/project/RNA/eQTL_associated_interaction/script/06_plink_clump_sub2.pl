#将../output/05_used_clump.txt 按照source，tissue和cis or trans进行分割，得../output/used_to_clump/${database}_${tissue}_cis.eSNP， ../output/used_to_clump/${database}_${tissue}_trans.eSNP
#../output/used_to_clump/${database}_${tissue}_unknown_cis_trans.eSNP 
#然后进行plink clunp,得../output/clump/${database}_${tissue}_cis.eSNP, ../output/clump/${database}_${tissue}_trans.eSNP, ../output/clump/${database}_${tissue}_unknown_cis_trans.eSNP
#进行PancanQTL的计算
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager; #多线程并行

my $fo1 = "../output/06_database_and_tissue_info_sub2.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Database\tTissue\n";

my %hash1;
my %hash2;
# my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
# my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource\tgene_type";
# print $O1 "$header\n";
#------------------------------------------------------ #取database and tissue
my @datasets = ("Blood_eQTL_browser","CAGE","eQTL_Catalog","eQTL_databases","eQTLGen","GTEx","PancanQTL","SEEQTL","ncRNA_eQTL");
my $Blood_eQTL_browser_tissue = "Blood";
my $CAGE_tissue = "Peripheral Blood";
my $eQTL_Catalog_tissue= "Blood";
my $eQTL_databases_tissue ="lymphoblastoid cell lines";
my $eQTLGen_tissue = "Blood";
my @GTEx_tissue;
my $f1 = "../data/GTEx/output/unique_tissue.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
while(<$I1>)
{
    chomp;
    push @GTEx_tissue,$_;
}

my @PancanQTL_tissue = ("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM");
my $SEEQTL_tissue = "NA";
my @ncRNA_eQTL_tissue =("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM");

my @datasets_sub1 = ("Blood_eQTL_browser","CAGE","eQTL_Catalog","eQTL_databases","eQTLGen","SEEQTL");
my @datasets_sub2 = ("GTEx","PancanQTL","ncRNA_eQTL");
# foreach my $dataset1(@datasets_sub1){
#     my $v = "${dataset1}_tissue";
#     push @{$hash1{$dataset1}},$v;
# }
# @{$hash1{GTEx}} = @GTEx_tissue;
@{$hash1{PancanQTL}} = @PancanQTL_tissue;
# @{$hash1{ncRNA_eQTL}} = @ncRNA_eQTL_tissue;

#------------------------------------------------------------------------------------------

foreach my $database(sort keys %hash1){
    my @tissues  = @{$hash1{$database}}; 
    foreach my $tissue(@tissues){
        # my $pid = $pm->start and next; #开始多线程
        print $O1 "$database\t$tissue\n";
        my $f2 = "../output/05_used_clump.txt";
        open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
        my $fo2 = "../output/used_to_clump/${database}_${tissue}_cis.eSNP";
        open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
        my $fo3 = "../output/used_to_clump/${database}_${tissue}_trans.eSNP";
        open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
        my $fo4 = "../output/used_to_clump/${database}_${tissue}_unknown_cis_trans.eSNP";
        open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";
        my $header ="SNP\tP";
        print $O2 "$header\n";
        print $O3 "$header\n";
        print $O4 "$header\n";

        while(<$I2>)
        {
            chomp;
            unless(/^Variant_id/){
                my @f =split/\t/;
                my $Rs_id = $f[0];
                unless($Rs_id =~/NA/){ #跳过rs id 为NA的行
                    my $Variant_Chr =$f[1];
                    $Variant_Chr =lc($Variant_Chr); #因为有些chr是以1表示，有的是以chr1表示，为了不进行条件判断而又都统一成1,所以先将chr替换掉，
                    $Variant_Chr =~s/chr//g;
                    my $chr ="chr${Variant_Chr}";
                    my $P_value =$f[2];
                    my $Q_value =$f[3];
                    my $Cis_or_Trans = $f[4];
                    $Cis_or_Trans = lc($Cis_or_Trans);
                    my $Tissue = $f[5];
                    my $Source =$f[6];
                    if ($Source =~/\b$database\b/ && $Tissue=~/\b$tissue\b/){
                        if ($database =~/SEEQTL/){  #SEEQTL 提供的是q value
                            if($Cis_or_Trans =~ /cis/){
                                print $O2 "$Rs_id\t$Q_value\n";
                            }
                            elsif($Cis_or_Trans =~ /trans/){
                                print $O3 "$Rs_id\t$Q_value\n";
                            }
                            else{ #不清楚trans和cis
                                print $O4 "$Rs_id\t$Q_value\n";
                            }
                        }
                        else{ #其余数据库提供的是P value
                            if($Cis_or_Trans =~ /cis/){
                                print $O2 "$Rs_id\t$P_value\n";
                            }
                            elsif($Cis_or_Trans =~ /trans/){
                                print $O3 "$Rs_id\t$P_value\n";
                            }
                            else{#不清楚trans和cis
                                print $O4 "$Rs_id\t$P_value\n";
                            }
                        }
                    }
                }
            }
        }
        close($I2);
        close($O2);
        close($O3);
        close($O4);
        my $command1 = "plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq --threads 20 --clump-r2 0.2 --clump-kb 500 --clump $fo2 --out ../output/clump/${database}_${tissue}_cis.eSNP";
        my $command2 = "plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq --threads 20 --clump-r2 0.2 --clump-kb 500 --clump $fo3 --out ../output/clump/${database}_${tissue}_trans.eSNP";
        my $command3 = "plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq --threads 20 --clump-r2 0.2 --clump-kb 500 --clump $fo4 --out ../output/clump/${database}_${tissue}_unknown_cis_trans.eSNP";
        
        my $line1s =readpipe("wc -l $fo2");
        my @c = split/\s+/,$line1s;
        my $line1 = $c[0];
        if ($line1 >1){  #如果cis文件有数据,进行$command1 
            system  $command1 ; #--threads 20 为多线程
        }
        #-----------------------       
        my $line2s =readpipe("wc -l $fo3");
        my @t = split/\s+/,$line2s;
        my $line2 = $t[0];
        if ($line2 >1){  #如果trnas文件有数据,进行$command2 
            system  $command2 ; #--threads 20 为多线程
        }       
#---------------------------------
        my $line3s =readpipe("wc -l $fo4");
        my @ss = split/\s+/,$line3s;
        my $line3 = $ss[0];
        if ($line3 >1){  #如果unknown_cis_trans文件有数据,进行$command3 
            system  $command3 ; #--threads 20 为多线程
        }        
        # print "$command1\n";
        # $pm->finish;  #多线程结束
    }
}
