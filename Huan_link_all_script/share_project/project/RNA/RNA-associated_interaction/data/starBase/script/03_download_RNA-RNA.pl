#download RNA-RNA
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @genetypes = ("mRNA","lncRNA","pseudogene","sncRNA","miRNA");

foreach my $geneType(@genetypes){
    my $command = "curl 'http://starbase.sysu.edu.cn/api/RNARNA/?assembly=hg19&geneType=${geneType}&RNA=all&interNum=1&expNum=1&cellType=all' > ../raw_data/ENCORI_hg19_${geneType}-RNA_RNA-RNA_all.txt";
    system $command;
}



# foreach my $geneType(@genetypes){
#     foreach my $program(@programs){
#         my $command="curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=${geneType}&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${program}&target=all&cellType=all' > ENCORI_hg19_CLIP-seq_all_miRNA_${geneType}_${program}_miRNA-Target.txt";
        
#     }
# }
