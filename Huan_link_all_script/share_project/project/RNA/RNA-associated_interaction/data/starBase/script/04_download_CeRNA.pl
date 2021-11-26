#download CeRNA
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @genetypes = ("mRNA","lncRNA","pseudogene");

foreach my $geneType(@genetypes){
    my $command = "curl 'http://starbase.sysu.edu.cn/api/ceRNA/?assembly=hg19&geneType=${geneType}&ceRNA=all&miRNAnum=2&family=all&pval=0.01&fdr=0.01&pancancerNum=0' > ../raw_data/ENCORI_hg19_ceRNA-${geneType}-network_all.txt";
    system $command;
}



# foreach my $geneType(@genetypes){
#     foreach my $program(@programs){
#         my $command="curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=${geneType}&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${program}&target=all&cellType=all' > ENCORI_hg19_CLIP-seq_all_miRNA_${geneType}_${program}_miRNA-Target.txt";
        
#     }
# }
