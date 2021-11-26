#用../output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt中的kinase提取
#../output/02_arrange_gene_effect_to_normal_type.txt中的kinase
#得kinase文件../output/04_kinase_Achilles.txt
#得dependent_probability 的kinase文件 ../output/04_kinase_Achilles_dependent.txt
#得非 dependent_probability 的kinase文件 ../output/04_kinase_Achilles_non_dependent.txt
#得非kinase 文件 ../output/04_non_kinase_Achilles.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../output/02_arrange_gene_effect_to_normal_type.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 = "../output/04_kinase_Achilles.txt"; #得kinase文件
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../output/04_kinase_Achilles_dependent.txt"; #得dependent_probability 的kinase文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 = "../output/04_kinase_Achilles_non_dependent.txt"; #非 dependent_probability 的kinase文件
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my $fo4 = "../output/04_non_kinase_Achilles.txt"; #得非kinase 文件
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";
my $fo5 = "../output/04_kinase_not_in_Achilles.txt"; ##不在Achilles中存在的kinase
open my $O5, '>', $fo5 or die "$0 : failed to open output file '$fo5' : $!\n";
my (%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7,%hash8,%hash9,%hash10,%hash11,%hash12,%hash13,%hash30);

my $header = "Gene\tSYMBOL\tENSG_id\tEntrez\tcell_lines(Broad_IDs)\tCERES";
print $O1 "$header\n";
print $O2 "$header\n";
print $O3 "$header\n";
print $O5 "Gene\tSYMBOL\tENSG_id\tEntrez\n";

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
        my $v ="$gene\t$SYMBOL\t$ensg";
        $hash1{$Entrez}=$v;
        # print "$ensg\n";
    }

}


while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    if (/^cell_lines/){
        print $O4 "$_\n";
    }
    else{
        my $cell_lines =$f[0];
        my $gene =$f[1];
        my $dependency_probability =$f[2];
        $gene =~ /(\(.*\))/; #捕获括号中的元素
        my $Entrez = $1;
        $Entrez =~ s/\(//g;
        $Entrez =~ s/\)//g;
        $hash2{$Entrez}=1;
        # print "$Entrez\t$gene\n";
        if (exists $hash1{$Entrez}){
            my $v = $hash1{$Entrez};
            my $output = "$v\t$Entrez\t$cell_lines\t$dependency_probability";
            print $O1 "$output\n";
            if ($dependency_probability =~/NA/){ #dependency_probability为NA
                print $O3 "$output\n";
            }
            else{
                if($dependency_probability < -0.6){ #dependency_probability >0
                    print $O2 "$output\n";
                }
                else{
                    print $O3 "$output\n"; #dependency_probability <0
                }
            }
        }
        else{
            print $O4 "$_\n";
        }
    }
}


foreach my $Entrez (sort keys %hash1){
    unless(exists $hash2{$Entrez}){  #不在Achilles中存在的kinase
        my $v = $hash1{$Entrez};
        print $O5 "$v\t$Entrez\n";
    }
}