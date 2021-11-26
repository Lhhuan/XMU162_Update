#download miRNA-Target
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @genetypes = ("mRNA", "lncRNA", "pseudogene", "circRNA", "sncRNA");
my @programs =("PITA,RNA22,miRmap","DIANA-microT,miRanda","PicTar,TargetScan");
foreach my $geneType(@genetypes){
    foreach my $program(@programs){
        my $command="curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=${geneType}&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${program}&target=all&cellType=all' > ../raw_data/ENCORI_hg19_CLIP-seq_all_miRNA-${geneType}_${program}_miRNA-Target.txt";
        system $command;
    }
}
