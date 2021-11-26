#为../output/01_all_kinds_QTL.txt添加genetype 得../output/02_annotate_gene_type.txt
#用"/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt"， "/home/huanhuan/ref_data/RNA_position/data/hg19_v19_miRNA_position.txt"，
#"/home/huanhuan/ref_data/RNA_position/data/hg19_mRNA_position.txt"和"/home/huanhuan/ref_data/RNA_position/data/gencode.v32lift37.long_noncoding_RNAs.gff3.gz"为../output/01_all_kinds_QTL.txt 注释gene type
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "/home/huanhuan/ref_data/RNA_position/data/hg19_v19_miRNA_position.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $f3 = "/home/huanhuan/ref_data/RNA_position/data/hg19_mRNA_position.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
my $f4 = "/home/huanhuan/ref_data/RNA_position/data/gencode.v32lift37.long_noncoding_RNAs.gff3.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
my $f5 = "../output/01_all_kinds_QTL.txt";
open my $I5, '<', $f5 or die "$0 : failed to open input file '$f5' : $!\n";
# open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件
my $fo1 = "../output/02_annotate_gene_type.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;
my %hash2;
print $O1 "unique_QTL_ID\tSNP_chr\tSNP_pos\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tPvalue\tSourceid\tQTL_type\tGene_type\n";


while(<$I1>)
{
    chomp;
    unless(/^Gene/){
        my @f =split/\t/;
        my $gene_name =$f[6];
        my $Gene_type = $f[7];
        $Gene_type =~ s/protein_coding/mRNA/g;
        $gene_name =uc($gene_name);
        # my $v = "$chr\t$start\t$end\t$Gene_type";
        $hash1{$gene_name}=$Gene_type;
    }
}


while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f =split/\t/;
        my $chr = $f[0];
        my $mi_type= $f[2];
        my $start =$f[3];
        my $end =$f[4];
        my $info =$f[-1];
        # my $v = "$mi_type\t$chr\t$start\t$end";
        # $
        my @ts = split/\;/,$info;
        foreach my $t(@ts){
            if ($t =~/Name=/){ 
                $t =~ s/Name=//g;
                $t = uc($t);
                unless(exists $hash1{$t}){
                    $hash1{$t}="miRNA";
                }
            }
        }
    }
}

while(<$I3>)
{
    chomp;
    unless(/^id/){
        my @f =split/\t/;
        my $chr = $f[2];
        my $start =$f[3];
        my $end =$f[4];
        my $gene_name =$f[1];
        $gene_name =uc($gene_name);
        # my $v = "$chr\t$start\t$end\tmRNA";
        my $v = "mRNA";
        unless(exists $hash1{$gene_name}){
            $hash1{$gene_name}=$v;
        }
    }
}


while(<$I4>)
{
    chomp;
    unless(/^#/){
        my @f =split/\s+/;
        my $chr = $f[0];
        my $gene_type =$f[2];
        my $start =$f[3];
        my $end =$f[4];
        my $info = $f[-1];
        if ($gene_type =~/gene/){
            my @ts = split/\;/,$info;
            foreach my $t(@ts){
                if ($t =~/gene_name/){
                    $t =~s/gene_name=//g;
                    $t =uc($t);
                    # my $v = "$chr\t$start\t$end\tlincRNA";
                    my $v = "lincRNA";
                    unless(exists $hash1{$t}){
                        $hash1{$t}=$v;
                    }
                }
            }
        }
    }
}

my $i =0;
while(<$I5>)
{
    chomp;
    unless(/^SNP_chr/){
        $i++;
        my @f =split/\t/;
        my $gene_name = $f[2];
        $gene_name =uc ($gene_name);
        if ( exists $hash1{$gene_name}){
            my $gene_type = $hash1{$gene_name};
            print $O1 "$i\t$_\t$gene_type\n";
        }
        else{ #没有匹配到gene type
            print $O1 "$i\t$_\tNULL\n";
        }
    }
}

close $O1;

foreach my $k (sort keys %hash1){
    my $v = $hash1{$k};
    print "$k\t$v\n";
}

    