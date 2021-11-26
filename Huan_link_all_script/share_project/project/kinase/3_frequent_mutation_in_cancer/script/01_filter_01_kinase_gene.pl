#从/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/mc3.v0.2.8.PUBLIC.maf.gz 中筛选出
#../../2_strong_gene-impact/output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt基因相关信息，得../output/01_TCGA_gene_in_01_kinase.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../../2_strong_gene-impact/output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/mc3.v0.2.8.PUBLIC.maf.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/01_TCGA_gene_in_01_kinase.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^query/){
        my $gene =$f[0];
        my $Entrez = $f[1];
        my $SYMBOL =$f[2];
        my $ensg = $f[3];
        $ensg =~ s/\s+//g;
        $ensg =~ s/list//g;
        $ensg =~ s/gene=//g;
        $ensg =~ s/,.*//g;
        $ensg =~ s/\(//g;
        $ensg =~ s/\)//g;
        $ensg =~ s/c//g;
        $ensg =~ s/"//g;
        # my $v ="$gene\t$SYMBOL\t$Entrez\t$ensg";
        # $hash1{$Entrez}=$v;
        # # print "$ensg\n";
        my $v = "$SYMBOL\t$Entrez\t$gene";
        $hash1{$ensg}=$v;
    }
}


while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    my $Hugo_Symbol =$f[0];
    my $NCBI_Build =$f[3];
    my $Chr = $f[4];
    my $start = $f[5];
    my $end= $f[6];
    my $Strand = $f[7];
    my $Variant_Classification =$f[8];
    my $Variant_Type =$f[9];
    my $ref = $f[10];
    my $tumor_allele1 = $f[11];
    my $tumor_allele2 = $f[12];
    my $dbSNP_RS =$f[13];
    my $Tumor_Sample_Barcode = $f[15];
    my $HGVSc =$f[34];
    my $HGVSp =$f[35];
    my $ensg =$f[47];
    my $output = "$Hugo_Symbol\t$NCBI_Build\t$Chr\t$start\t$end\t$Strand\t$Variant_Classification\t$Variant_Type\t$ref\t$tumor_allele1";
    $output = "$output\t$tumor_allele2\t$dbSNP_RS\t$Tumor_Sample_Barcode\t$HGVSc\t$HGVSp\t$ensg";
    if (/^Hugo_Symbol/){
        print $O1 "$output\tSYMBOL\tEntrez\tgene\n";
    }
    else{
        if (exists $hash1{$ensg}){
            my $v = $hash1{$ensg};
            print $O1 "$output\t$v\n";
        }
    }
}
