#download data 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @array = ("ACC","BLCA","BRCA","CESC","COAD","DLBC","ESCA","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM");
for  my $cancer (@array){ 
# my $link = "http://ibi.hzau.edu.cn/ncRNA-eQTL/miRNA/download_data/gwas/${cancer}_GWAS_eQTLs.txt";
my $link = "http://ibi.hzau.edu.cn/ncRNA-eQTL/download_data/survival/${cancer}_Survival_eQTLs.txt";
my $commod = "wget ${link}\n";
system "$commod";
}