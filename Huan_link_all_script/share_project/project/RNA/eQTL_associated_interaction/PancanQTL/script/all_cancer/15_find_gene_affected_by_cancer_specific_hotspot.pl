#找 hotspot对应的eGene $output_dir/cancer_${cancer}_${tissue}_specific.bed.gz 和"$output_dir/qtl_gene.bed.gz" 得 "$output_dir/cancer_${cancer}_${tissue}_specific_gene.bed.gz"
#得汇总文件"../../output/cancer_total/15_cancer_specfic_hotspot_gene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;

my $f1 = "../../data/pancanQTL_gtex_eQTL.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^Study/){
       my @f= split/\t/;
       my $cancer=$f[0];
       $cancer =~ s/\s+//g;
       my $GTEx_tissue=$f[-1];
       my @ts=split/;/,$GTEx_tissue;
    #    print "$cancer\n";
       foreach my $t(@ts){
           push @{$hash1{$cancer}},$t;
       }
    }
}    

# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","UCEC","UCS"); 
foreach my $cancer(@cancers){
    # print "$cancer\n";
    if(exists $hash1{$cancer}){
        # print "$cancer\n";
        my @tissues =@{$hash1{$cancer}};
        foreach my $tissue(@tissues){
            my $output_dir = "../../output/cancer_total/specific/pure/${cancer}";
            my $cancer_sp = "$output_dir/cancer_${cancer}_${tissue}_specific.bed.gz";
            my $gene = "$output_dir/qtl_gene.bed.gz";
            my $cancer_gene_sp ="$output_dir/cancer_${cancer}_${tissue}_specific_gene.bed.gz";
            my $command1 = "bedtools intersect -a $cancer_sp -b $gene -wa -wb |gzip >$cancer_gene_sp";
            # print "$command1\n";
            system $command1;
            print "$cancer\t$tissue\n";

        }
    }
}
#THCA 不明原因报错，故解压后操作
system "bash 15_follow_bedtools.sh";

my $fo1 = "../../output/cancer_total/15_cancer_specfic_hotspot_gene.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "h_chr\th_start\th_end\teqtl_chr\teqtl_start\teqtl_end\tegene\tpvalue\tcancer\ttissue\n";

my @cancer2s =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
foreach my $cancer(@cancer2s){
    if (exists $hash1{$cancer}){
        my @tissues =@{$hash1{$cancer}};
        foreach my $tissue(@tissues){
            my $output_dir = "../../output/cancer_total/specific/pure/${cancer}";
            my $cancer_gene_sp ="$output_dir/cancer_${cancer}_${tissue}_specific_gene.bed.gz";
            my $f2 = $cancer_gene_sp;
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件 
            while(<$I2>)
            {
                chomp;
            print $O1 "$_\t$tissue\n";
        
            }    
        }
    }
}


