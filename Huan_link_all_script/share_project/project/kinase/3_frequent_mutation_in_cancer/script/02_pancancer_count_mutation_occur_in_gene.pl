#pancancer 统计../output/01_TCGA_gene_in_01_kinase.txt.gz中 每个gene的突变情况，
#得../output/02_pancancer_count_mutation_occur_in_gene.txt
#并同时根据 /share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/TCGA-Clinical-Data-Resource.csv
#patients对应的具体cancer_type信息，得../output/02_TCGA_gene_in_01_kinase_cancer.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/TCGA-Clinical-Data-Resource.csv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../output/01_TCGA_gene_in_01_kinase.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/02_TCGA_gene_in_01_kinase_cancer.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/02_pancancer_count_mutation_occur_in_gene.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my (%hash1,%hash2,%hash3);
# my %hash2;
print $O2 "SYMBOL\tEntrez\tgene\tnumber\n";
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $bcr_patient_barcode =$f[1];
    my $type =$f[2];
    my $histological_type =$f[8];
    unless($bcr_patient_barcode =~ /bcr_patient_barcode/){
        my $v = "$type\t$histological_type";
        $hash1{$bcr_patient_barcode}=$v;
        # print "$bcr_patient_barcode\n";
    }
}


while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    if (/^Hugo_Symbol/){
        print $O1 "$_\tcancer_type\thistological_type\n";
    }
    else{
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
        my @t =split/\-/,$Tumor_Sample_Barcode;
        my $new_Tumor_Sample_Barcode  = join("-",@t[0..2]);
        print "$new_Tumor_Sample_Barcode\n";
        if (exists $hash1{$new_Tumor_Sample_Barcode}){
            my $v= $hash1{$new_Tumor_Sample_Barcode};
            print $O1 "$_\t$v\n";
        }
        my $variants = "$ensg\t$HGVSc\t$new_Tumor_Sample_Barcode";
        my $k2 = "$SYMBOL\t$Entrez\t$gene";
        push @{$hash2{$k2}},$variants;
    }
}

foreach my $k2(sort keys %hash2){
    my @v2s = @{$hash2{$k2}};
    @v2s = grep { ++$hash3{$_} < 2 } @v2s;
    my $number = @v2s;
    print $O2 "$k2\t$number\n";
    # my $ou = join("\n",@v2s);
    # print "$ou\t$k2\n";
}