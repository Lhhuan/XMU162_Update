#transform goa_human.gaf.gz 的symbol to ensg得 /home/huanhuan/project/eQTL_Catalogue/script/prediction/Edgeprediction/output/train/03_go_edgepredict.csv
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "goa_human.gaf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "02_symbol_ENSG.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $fo1 = "/home/huanhuan/project/eQTL_Catalogue/script/prediction/Edgeprediction/output/train/03_go_edgepredict.csv"; #
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;

# my $fo2 = "../output/train/01_reactome.csv"; #
# open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";


print $O1 "Source node name,Source node type,Relationship type,Target node type,Target node name\n";

my (%hash1,%hash2);
while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    unless(/^query/){
        my $gene =$f[0];
        my $ensg = $f[-4];
        # my $Entrez = $f[-2];
        # my $SYMBOL =$f[-1];
        $ensg =~ s/\s+//g;
        $ensg =~ s/list//g;
        $ensg =~ s/gene=//g;
        $ensg =~ s/,.*//g;
        $ensg =~ s/\(//g;
        $ensg =~ s/\)//g;
        $ensg =~ s/c//g;
        $ensg =~ s/"//g;
        # my $v ="$Entrez\t$SYMBOL\t$ensg";
        # push @{$hash1{$gene}},$v;
        push @{$hash2{$gene}},$ensg;
    }
}


while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^\!/){
        my $gene=$f[2];
        my $term =$f[4];
        if(exists $hash2{$gene}){
            my @v2s=@{$hash2{$gene}};
            my $v2=$v2s[0];
            print $O1 "$v2,gene,go_annotation,GOterm,$term\n"

        }
    }
}