#得"../../output/cancer_total/08_number_of_cancer_all_hotspot.txt.gz"，"../../output/cancer_total/08_number_of_cancer_specific_hotspot.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $cancer = "Lung";
my $j = 18;

# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 

my @cancers =("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM");
 
my $fo1 = "../../output/cancer_total/trans/05_number_of_cancer_all_hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;


print $O1 "cancer\tNumber\n";

 
foreach my $cancer(@cancers){
    
    my $f1 = "../../output/${cancer}/Trans_eQTL/hotspot_trans_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    if(-e $f1 ){
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        my %hash1;
        while(<$I1>)
        {
            $hash1{$_}=1;
        }
        my $hash_length = keys %hash1;
        print $O1 "$cancer\t$hash_length\n";
    }
    
}


sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
