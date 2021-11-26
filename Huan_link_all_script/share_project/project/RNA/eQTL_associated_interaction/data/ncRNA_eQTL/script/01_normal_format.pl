#将../raw_data/Cis-eQTLs.txt 和../raw_data/Trans-eQTLs.txt 进行初步normal, ../output/01_normal_format.txt.gz
#并得unique的gene POS的bed文件../output/01_unique_gene_position_hg38.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "../output/01_normal_format.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/01_unique_gene_position_hg38.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $header = "Variant_id\tRs_id\tChr\tPos\tRef\tAlt\tMaf\tGene_Version\tGene_id\tGene_name\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue";
print $O1 "$header\n";
print $O2 "pos\tstart\tend\tgene_id\n";
my @types = ("Cis","Trans");
my(%hash1,%hash2,%hash3,%hash4);
foreach my $type(@types){
    my $f1 = "../raw_data/${type}-eQTLs.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
    # print "$f1\n";
    while(<$I1>){
        chomp;
        my @f =split/\t/;
        unless(/^cancer_type/){
            my $Variant_id = $f[1];
            my $Rs_id =$f[1];
            my $SNP_position =$f[9];
            my @c =split/\:/,$SNP_position;
            my $Chr=$c[0];
            # my $Pos= $c[1];
            my $Pos= "NA"; #此处不再记录hg38的pos,等后面一起用hg19的annotation进行填充
            my $Alleles = $f[10];
            my @t =split/\//,$Alleles;
            my $Ref= $t[0];
            my $Alt= $t[1];
            # print "$Ref\t$Alt\n";
            my $Maf= "NA";
            my $Gene_Version="Hg38";
            my $Gene_id= $f[2];
            my $Gene_name= $f[3] ;
            my $Slope= "NA";
            my $Slope_se= "NA";
            my $P_value=$f[4];
            my $Q_value= "NA";
            my $Beta= $f[7];
            my $Tissue=$f[0] ;
            my $Effect_size = "NA";
            my $Zscore= "NA";
            my $Cis_or_Trans=$type ;
            my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effect_size\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue";
            print $O1 "$output\n";
            my $Gene_position =$f[-5];
            $Gene_position =~ s/~/:/g;
            my @s = split/\:/,$Gene_position;
            my $gene_chr = $s[0];
            my $gene_start = $s[1];
            my $gene_end = $s[2];
            my $gene_pos_info = "$gene_chr\t$gene_start\t$gene_end\t$Gene_id";
            unless (exists $hash3{$gene_pos_info}){
                $hash3{$gene_pos_info} =1;
                # print "$gene_pos_info\n";
                print $O2 "$gene_pos_info\n";
            }
        }
    }
    close $I1;
}
close $O1;
