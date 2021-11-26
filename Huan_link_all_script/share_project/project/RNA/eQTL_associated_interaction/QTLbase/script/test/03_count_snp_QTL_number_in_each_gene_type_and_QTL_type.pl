#count ../output/02_annotate_gene_type.txt.gz 各类QTL的数量,snp的数量及基因的数量得文件../output/03_QTL_type_number.txt  在特定QTL种类下的QTL-gene type interaction文件得../output/03_QTL_type_gene_type_number.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../output/02_annotate_gene_type.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/03_QTL_type_number.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "../output/03_QTL_type_gene_type_number.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

my %hash1;
my %hash2;
my %hash3;
# print $O1 "unique_QTL_ID\tSNP_chr\tSNP_pos\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tPvalue\tSourceid\tQTL_type\tGene_type\n";
print $O1 "QTL_type\tQTL_number\tSNP_number\tGene_number\n";
print $O2 "QTL_type\tinteraction_Gene_type\tSNP_number\n";

while(<$I1>)
{
    chomp;
    unless(/^unique_QTL_ID/){
        my @f =split/\t/;
        my $unique_QTL_ID = $f[0];
        my $SNP_chr =$f[1];
        my $SNP_pos =$f[2];
        my $Mapped_gene =$f[3];
        my $QTL_type =$f[-2];
        my $Gene_type =$f[-3];
        my $v1 = "$unique_QTL_ID\t$SNP_chr\t$SNP_pos\t$Mapped_gene";
        push @{$hash1{$QTL_type}},$v1; #各种type的QTL;
        my $k2 ="$QTL_type\t$Gene_type";
        my $v2 = "$Mapped_gene\t$unique_QTL_ID";
        push @{$hash2{$k2}},$v2;  #intraction type in QTL type
        my $k3 ="$Gene_type";
        my $v3 = "$Mapped_gene\t$unique_QTL_ID";        
        push @{$hash3{$k3}},$v3;  #all intraction type
    }
}
#-------------------------------------------------------count snp,gene in each QTL type
my $all_QTL_number =0;
my @all_snp;
my @all_gene;
foreach my $QTL_type(sort keys %hash1){
    my @v1s = @{$hash1{$QTL_type}};
    my $qtl_number = @v1s;
    $all_QTL_number = $all_QTL_number +$qtl_number; #count all QTL
    my @snps;
    my @genes;
    foreach my $v1(@v1s){
        my @f =split/\t/,$v1;
        my $Snp_position = join("\t",@f[1,2]);
        push @snps,$Snp_position; #  qtl specific
        push @all_snp,$Snp_position; #all
        my $gene= $f[-1];
        push @genes, $gene; #  qtl specific
        push @all_gene,$gene; #all
    }
    my %hash4;
    my %hash5;
    @snps = grep { ++$hash4{$_} < 2 } @snps;
    @genes = grep {++$hash5{$_}<2}@genes;
    my $snp_number = @snps;
    my $genes_number = @genes;
    print $O1 "$QTL_type\t$qtl_number\t$snp_number\t$genes_number\n";
}
my %hash6;
my %hash7;
@all_snp = grep { ++$hash6{$_} < 2 } @all_snp;
my $all_snp_number = @all_snp;
@all_gene = grep { ++$hash7{$_} < 2 } @all_gene;
my $all_gene_number = @all_gene;
print $O1 "Total\t$all_QTL_number\t$all_snp_number\t$all_gene_number\n";
#-----------------------------------------------------------------------------------------

#------------------------------------------------QTL type and gene type 
foreach my $k2 (sort keys %hash2){
    my @v2s = @{$hash2{$k2}};
    my $number = @v2s;
    print $O2 "$k2\t$number\n";
}

foreach my $k3(sort keys %hash3){
    my @v3s = @{$hash3{$k3}};
    my $number = @v3s;
    print $O2 "Total\t$k3\t$number\n";
}