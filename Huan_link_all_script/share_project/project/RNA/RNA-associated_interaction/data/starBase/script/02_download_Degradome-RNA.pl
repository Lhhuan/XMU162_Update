#download Degradome-RNA
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @genetypes = ("mRNA","ncRNA");

foreach my $geneType(@genetypes){
    my $command = "curl 'http://starbase.sysu.edu.cn/api/degradomeRNA/?assembly=hg19&geneType=${geneType}&miRNA=all&degraExpNum=1&target=all&cellType=all' > ../raw_data/ENCORI_hg19_degradome-seq_all_Degradome-${geneType}_Degradome-RNA.txt";
    system $command;
}



# foreach my $geneType(@genetypes){
#     foreach my $program(@programs){
#         my $command="curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=${geneType}&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${program}&target=all&cellType=all' > ENCORI_hg19_CLIP-seq_all_miRNA_${geneType}_${program}_miRNA-Target.txt";
        
#     }
# }
