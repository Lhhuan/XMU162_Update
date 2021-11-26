#在../output/05_merge_essential_gene_cell_line_cancer.txt中筛选出有TCGA_code的essential gene
#得../output/06_tcga_code_essential_gene_annotation.txt 和 非tcga code 的文件../output/06_NON-tcga_code_essential_gene_annotation.txt
#得每个tcga code 下的essential_gene文件，../output/06_tcga_code_essential_gene.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/05_merge_essential_gene_cell_line_cancer.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../output/06_NON-tcga_code_essential_gene_annotation.txt"; #得kinase cancer 文件文件
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../output/06_tcga_code_essential_gene_annotation.txt"; #得kinase cancer 文件文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $fo3= "../output/06_tcga_code_essential_gene.txt"; #得kinase cancer 文件文件
open my $O3, '>', $fo3 or die "$0 : failed to open output file  '$fo3' : $!\n";
my $fo4= "../output/06_pancancer_essential_gene.txt"; #得pancancer的kinase 文件
open my $O4, '>', $fo4 or die "$0 : failed to open output file  '$fo4' : $!\n";
my %hash1;
my %hash3;

my $header = "Gene\tSYMBOL\tENSG_id\tEntrez\tcell_lines(Broad_IDs)\tCERES\ttcga_code\tdisease\tPATHOLOGIST_ANNOTATION\tsource";
print $O1 "$header\n";
print $O2 "$header\n";
print $O3 "TCGA_code\tSYMBOL\n";
print $O4 "gene\tSYMBOL\tENSG_id\tEntrez\n";

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^Gene/){
        my $gene =$f[0];
        my $SYMBOL =$f[1];
        my $ENSG_id =$f[2];
        my $Entrez =$f[3];
        my $cell_lines =$f[4];
        my $CERES = $f[5];
        my $tcga_code =$f[6];
        my $disease = $f[7];
        my $PATHOLOGIST_ANNOTATION =$f[8];
        my $source = $f[9];
        my $output_gene = "$gene\t$SYMBOL\t$ENSG_id\t$Entrez";
        unless(exists $hash3{$output_gene}){
            $hash3{$output_gene} =1;
            print $O4 "$output_gene\n";
        }
        if ($tcga_code =~/\bNA\b/){
            print $O1 "$_\n";
        }
        elsif($tcga_code =~ /UNABLE/){
            print $O1 "$_\n";
        }
        else{
            print $O2 "$_\n";
            push @{$hash1{$tcga_code}},$SYMBOL;
        }
    }

}

foreach my $tcga_code (sort keys %hash1){
    my @genes = @{$hash1{$tcga_code}};
    my %hash2;
    @genes =grep {++$hash2{$_}<2} @genes;
    my $output_gene = join(";",@genes);
    print $O3 "$tcga_code\t$output_gene\n";
}