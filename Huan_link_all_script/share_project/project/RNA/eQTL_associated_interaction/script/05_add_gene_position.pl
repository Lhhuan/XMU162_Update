#用"/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt"， "/home/huanhuan/ref_data/RNA_position/data/hg19_v19_miRNA_position.txt"，
#"/home/huanhuan/ref_data/RNA_position/data/hg19_mRNA_position.txt"和"/home/huanhuan/ref_data/RNA_position/data/gencode.v32lift37.long_noncoding_RNAs.gff3.gz"为../output/04_add_gene_ENSG_and_Entrez.txt.gz 寻找gene的
#位置信息，得../output/05_add_gene_position.txt.gz
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
my $f5 = "../output/04_add_gene_ENSG_and_Entrez.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件
my $fo1 = "../output/05_add_gene_position.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;
my %hash2;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource\tgene_type";
print $O1 "$header\n";


while(<$I1>)
{
    chomp;
    unless(/^Gene/){
        my @f =split/\t/;
        my $chr = $f[3];
        my $start =$f[4];
        my $end =$f[5];
        my $gene_name =$f[6];
        my $Gene_type = $f[7];
        $Gene_type =~ s/protein_coding/mRNA/g;
        $gene_name =uc($gene_name);
        my $v = "$chr\t$start\t$end\t$Gene_type";
        $hash1{$gene_name}=$v;
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
        my $v = "$mi_type\t$chr\t$start\t$end";
        my @ts = split/\;/,$info;
        foreach my $t(@ts){
            if ($t =~/Name=/){ 
                $t =~ s/Name=//g;
                $t = uc($t);
                push @{$hash2{$t}},$v;

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
        my $v = "$chr\t$start\t$end\tmRNA";
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
                    my $v = "$chr\t$start\t$end\tlincRNA";
                    unless(exists $hash1{$t}){
                        $hash1{$t}=$v;
                    }
                }
            }
        }
    }
}


while(<$I5>)
{
    chomp;
    unless(/^Variant_id/){
        my @f =split/\t/;
        my $Gene_name = $f[11];
        my $EntrezGeneID =$f[12];
        my $chr = $f[13];
        my $start = $f[14];
        my $end =$f[15];
        my $Variant_info = join("\t",@f[0..10]);
        my $other_info = join("\t",@f[16..$#f]);
        my $gene_name = $Gene_name;
        $gene_name =uc ($gene_name);
        my $pos = "$chr\t$start\t$end";
        if ($pos =~/NA/){#chr start end 有为NA的项
            if(exists $hash1{$gene_name}){
                my $vs= $hash1{$gene_name};
                my @t = split /\t/,$vs;
                my $gene_type = $t[-1];
                my $v = join ("\t",@t[0..2]);
                my $output = "$Variant_info\t$Gene_name\t$EntrezGeneID\t$v\t$other_info\t$gene_type"; #hash1的v中记录了gene type
                print $O1 "$output\n";
            }
            elsif(exists $hash2{$gene_name}){ #hash2是miRNA
                my @vs = @{$hash2{$gene_name}};
                my %hash3;
                @vs = grep { ++$hash3{ $_ } < 2; } @vs;
                my $number = @vs;
                if ($number <2){ #只有一个v;
                    my $v = $vs[0];
                    my @t = split/\t/,$v;
                    my $pos = join ("\t",@t[1..3]);
                    my $output = "$Variant_info\t$Gene_name\t$EntrezGeneID\t$pos\t$other_info\tmiRNA"; 
                    print $O1 "$output\n";
                }
                else{#两个v
                    foreach my $v(@vs){
                        unless ($v =~/miRNA_primary_transcript/){#本文件中一共有8个miRNA,这8个具有多个v的mirna中分别对应一个miRNA_primary_transcript和一个miRNA,
                            my @t = split/\t/,$v;
                            my $pos = join ("\t",@t[1..3]);
                            my $output = "$Variant_info\t$Gene_name\t$EntrezGeneID\t$pos\t$other_info\tmiRNA";
                            print $O1 "$output\n";                           
                        }
                    }
                }
            }
            else{ #gene name 在hash1和hash2中都不存在，按照原文输出,gene type记录为NA
                print $O1 "$_\tNA\n";
            }
        }
        else{ #chr start end 均不为NA, 按照原文输出
            if(exists $hash1{$gene_name}){
                my $vs= $hash1{$gene_name};
                my @t = split /\t/,$vs;
                my $gene_type = $t[-1];
                print $O1 "$_\t$gene_type\n";
            }
            elsif(exists $hash2{$gene_name}){ #hash2是miRNA
                print $O1 "$_\tmiRNA\n";
            }
            else{ #gene name 在hash1和hash2中都不存在，按照原文输出,gene type记录为NA
                print $O1 "$_\tNA\n";
            }
        }
    }
}

close $O1;

    