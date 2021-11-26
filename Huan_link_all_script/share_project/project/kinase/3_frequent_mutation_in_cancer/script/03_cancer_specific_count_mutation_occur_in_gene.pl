#cancer specific 统计../output/02_TCGA_gene_in_01_kinase_cancer.txt中 每个gene的突变情况，
#得../output/03_cancer_specific_count_mutation_occur_in_gene.txt 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/02_TCGA_gene_in_01_kinase_cancer.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../output/03_cancer_specific_count_mutation_occur_in_gene.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my (%hash1,%hash2,%hash3);
# my %hash2;
print $O1 "SYMBOL\tEntrez\tgene\tcancer_type\tnumber\n";


while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless (/^Hugo_Symbol/){
        my $Hugo_Symbol =$f[0];
        my $NCBI_Build =$f[1];
        my $Chr =$f[2];
        my $start = $f[3];
        my $end =$f[4];
        my $Strand =$f[5];
        my $Variant_Classification =$f[6];
        my $Variant_Type = $f[7];
        my $ref =$f[8];
        my $tumor_allele1 =$f[9];
        my $tumor_allele2 =$f[10];
        my $dbSNP_RS =$f[11];
        my $Tumor_Sample_Barcode =$f[12];
        my $HGVSc =$f[13];
        my $HGVSp =$f[14];
        my $ensg =$f[15];
        my $SYMBOL= $f[16];
        my $Entrez= $f[17];
        my $gene = $f[18];
        my $cancer_type =$f[19];
        my $histological_type =$f[20];
        my @t =split/\-/,$Tumor_Sample_Barcode;
        my $new_Tumor_Sample_Barcode  = join("-",@t[0..2]);
        my $variants = "$ensg\t$HGVSc\t$new_Tumor_Sample_Barcode";
        my $k2 = "$SYMBOL\t$Entrez\t$gene\t$cancer_type";
        push @{$hash2{$k2}},$variants;
    }
}

foreach my $k2(sort keys %hash2){
    my @v2s = @{$hash2{$k2}};
    @v2s = grep { ++$hash3{$_} < 2 } @v2s;
    my $number = @v2s;
    print $O1 "$k2\t$number\n";
    # my $ou = join("\n",@v2s);
    # print "$ou\t$k2\n";
}