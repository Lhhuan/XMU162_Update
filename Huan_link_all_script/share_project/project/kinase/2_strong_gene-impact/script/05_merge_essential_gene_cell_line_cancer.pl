#用../data/DepMap/Cell_lines_annotations_20181226.txt 和../data/DepMap/sample_info.csv为../output/04_kinase_Achilles_dependent.txt
#寻找cell line的信息。首先用../data/DepMap/Cell_lines_annotations_20181226.txt寻找，然后用../output/04_kinase_Achilles_dependent.txt寻找
#得../output/05_merge_essential_gene_cell_line_cancer.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../data/DepMap/Cell_lines_annotations_20181226.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../data/DepMap/sample_info.csv";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $f3 = "../output/04_kinase_Achilles_dependent.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
my $fo1 = "../output/05_merge_essential_gene_cell_line_cancer.txt"; #得kinase cancer 文件文件
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my (%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7,%hash8,%hash9,%hash10,%hash11,%hash12,%hash13,%hash30);

my $header = "Gene\tSYMBOL\tENSG_id\tEntrez\tcell_lines(Broad_IDs)\tCERES\ttcga_code\tdisease\tPATHOLOGIST_ANNOTATION\tsource";
print $O1 "$header\n";

while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^CCLE_ID/){
        my $depMapID =$f[1];
        my $Disease =$f[18];
        my $PATHOLOGIST_ANNOTATION =$f[-3];
        my $tcga_code =$f[-1];
        my $v = "$tcga_code\t$Disease\t$PATHOLOGIST_ANNOTATION";
        $hash1{$depMapID}=$v;
        # print "$PATHOLOGIST_ANNOTATION\n";
    }

}
while(<$I2>)
{
    chomp;
    my @f= split/\,/;
    unless(/^DepMap_ID/){
        my $depMapID =$f[0];
        my $Disease =$f[5];
        my $PATHOLOGIST_ANNOTATION ="NA";
        my $tcga_code ="NA";
        my $v = "$tcga_code\t$Disease\t$PATHOLOGIST_ANNOTATION";
        $hash2{$depMapID}=$v;
    }

}

while(<$I3>)
{
    chomp;
    my @f= split/\t/;
    unless(/^Gene/){
        my $Gene =$f[0];
        my $SYMBOL =$f[1];
        my $ensg = $f[2];
        my $Entrez = $f[3];
        my $cell_lines = $f[4];
        my $CERES =$f[5];
        if (exists $hash1{$cell_lines}){ #
            my $v = $hash1{$cell_lines};
            print $O1 "$_\t$v\tCell_lines_annotations\n"; #来自Cell_lines_annotations的annotation
        }
        elsif(exists $hash2{$cell_lines}){
            my $v = $hash2{$cell_lines};
            print $O1 "$_\t$v\tsample_info\n"; #来自sample_info的annotation
        }
        else{
            print "$cell_lines\n";
        }
    }
}


