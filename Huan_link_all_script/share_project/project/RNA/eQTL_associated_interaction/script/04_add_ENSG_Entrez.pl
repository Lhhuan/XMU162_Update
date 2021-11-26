#用"~/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt" 为../output/03_add_rsid.txt.gz添加 ENSG和Entrez,得../output/04_add_gene_ENSG_and_Entrez.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/huanhuan/ref_data/Ensembl/output/Hg19_reference_geneome_data_from_ensembl_BioMart_entrez.txt";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../output/03_add_rsid.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/04_add_gene_ENSG_and_Entrez.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;
my $header = "Variant_id\tRs_id\tVariant_Chr\tVariant_Pos\tVariant_Ref\tVariant_Alt\tEffect_allele\tOther_allele\tMaf\tGene_Version\tGene_id\tGene_name\tEntrezGeneID\tGene_Chr\tGene_Start\tGene_End\tSlope\tSlope_se\tZscore\tEffct_sise\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue\tSource";
print $O1 "$header\n";


while(<$I1>)
{
    chomp;
    unless(/^Gene/){
        my @f =split/\t/;
        my $ENSG = $f[0];
        my $Gene_name =$f[6];
        my $Entrezgene =$f[-2];
        $Gene_name =uc($Gene_name);
        my $v = "$ENSG\t$Entrezgene";
        $hash1{$Gene_name}=$v;
    }
}


while(<$I2>)
{
    chomp;
    unless(/^Variant_id/){
        my @f =split/\t/;
        my $Gene_ensg_id = $f[10];
        my $original_Gene_name = $f[11];
        my $ddd= "///";
        $original_Gene_name =~ s/$ddd/,/g; #把///替换为,
        my $EntrezGeneID =$f[12];
        my $Variant_info = join("\t",@f[0..9]);
        my $other_info = join("\t",@f[13..$#f]);
        my @Gene_names = split/\,/,$original_Gene_name; #原纪录中有些gene name 以,分割，如：AC000032.7,GSTM1,GSTM2
        foreach my $Gene_name(@Gene_names){
            my $gene_name = $Gene_name; #为了保证gene name 按照原样子输出
            $gene_name =uc($gene_name);
            my $all_gene_id ="$Gene_ensg_id\t$EntrezGeneID";
            if ($all_gene_id =~/NA/){ #ensg和Entrez中有NA
                if(exists $hash1{$gene_name}){
                    my $v = $hash1{$gene_name};
                    my @t = split/\t/,$v;
                    my $new_ensg = $t[0];
                    my $new_Entrez =$t[1];
                    if($Gene_ensg_id =~ /NA/ && $EntrezGeneID =~/NA/){ #ensg和Entrez都是NA，用hash1中的ensg和Entrez
                        my $output ="$Variant_info\t$new_ensg\t$Gene_name\t$new_Entrez\t$other_info"; #输出的gene name是没有进行大写变换的
                        print $O1 "$output\n";
                    }
                    elsif($Gene_ensg_id !~ /NA/ && $EntrezGeneID =~/NA/){#ensg不为NA,Entrez是NA，只用hash1中的Entrez
                        my $output ="$Variant_info\t$Gene_ensg_id\t$Gene_name\t$new_Entrez\t$other_info";
                        print $O1 "$output\n";
                    }
                    else{ ##ensg为NA,Entrez不是NA，只用hash1中的ensg
                        my $output ="$Variant_info\t$new_ensg\t$Gene_name\t$EntrezGeneID\t$other_info";
                        print $O1 "$output\n";
                    }
                }
                else{ #gene name 不在hash1,输出原文中的Entrez,ensg
                    my $output = "$Variant_info\t$Gene_ensg_id\t$Gene_name\t$EntrezGeneID\t$other_info";
                    print $O1 "$output\n";
                }
            }
            else{ #ensg和Entrez都不是NA,(gene name用，分割后)原文输出
                my $output = "$Variant_info\t$Gene_ensg_id\t$Gene_name\t$EntrezGeneID\t$other_info";
                print $O1 "$output\n";
            }
        }
    }
}

close $O1;

    