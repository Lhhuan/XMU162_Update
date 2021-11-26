#download RBP-Target
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @genetypes = ("mRNA","lncRNA","pseudogene","circRNA","sncRNA");

foreach my $geneType(@genetypes){
    my $command = "curl 'http://starbase.sysu.edu.cn/api/RBPTarget/?assembly=hg19&geneType=${geneType}&RBP=all&clipExpNum=5&pancancerNum=0&target=all&cellType=all' > ../raw_data/ENCORI_hg19_RBP-${geneType}_RBP-Target.txt";
    system $command;
}



# foreach my $geneType(@genetypes){
#     foreach my $program(@programs){
#         my $command="curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=${geneType}&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${program}&target=all&cellType=all' > ENCORI_hg19_CLIP-seq_all_miRNA_${geneType}_${program}_miRNA-Target.txt";
        
#     }
# }
