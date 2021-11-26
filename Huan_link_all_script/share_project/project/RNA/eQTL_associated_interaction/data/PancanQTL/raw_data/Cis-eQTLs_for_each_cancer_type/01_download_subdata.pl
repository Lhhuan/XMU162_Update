#download data 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @array = ("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM");
for  my $cancer (@array){ 
my $link = "http://bioinfo.life.hust.edu.cn/PancanQTL/static/download/${cancer}_tumor.cis_eQTL.xls";
my $commod = "wget ${link}\n";
system "$commod";
}