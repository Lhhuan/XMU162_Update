##将../raw_data/MRCA_5percentFDR.csv.gz 和../raw_data/MRCE_5percentFDR.csv.gz 进行初步normal, 得gene name 和部分RSid从对应的v8.egenes.txt.gz中获得，得文件../output/01_normal_format.txt.gz
#并得unique的varint POS的vcf文件../output/01_unique_variant-id_position_hg18.vcf.gz
#因为不清楚../raw_data/meta_5percentFDR.csv.gz 所以不将该文件纳入
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "../output/01_normal_format.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../output/01_unique_variant-id_position_hg18.vcf.gz";
open my $O2, "| gzip >$fo2" or die $!;
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $header = "Variant_id\tRs_id\tChr\tPos\tRef\tAlt\tMaf\tGene_Version\tGene_id\tGene_name\tSlope\tSlope_se\tZscore\tEffect_size\tP_value\tQ_value\tBeta\tCis_or_Trans\tTissue";
print $O1 "$header\n";
print $O2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
my(%hash1,%hash2,%hash3,%hash4);
my $f1 = "../raw_data/MRCE_5percentFDR.csv.gz";
my $f2 = "../raw_data/MRCA_5percentFDR.csv.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
# print "$f1\n";
while(<$I1>){
    chomp;
    my @f =split/\,/;
    unless(/^gene/){
        my $Variant_id = $f[2];
        my $Rs_id = "NA";
        my $Chr=$f[3];
        my $Pos= $f[4];
        $Pos =~ s/\.//g;
        $Pos =~ s/\.//g;
        my $Ref= $f[5];
        my $Alt= $f[6];
        # print "$Ref\t$Alt\n";
        my $Maf= "NA";
        my $Gene_Version= "Hg18";
        my $Gene_id= "NA";
        my $Slope="NA";
        my $Slope_se= $f[11];
        my $P_value="NA";
        my $Q_value= "NA";
        my $Beta="NA";
        my $Tissue="lymphoblastoid cell lines";
        my $Gene_name= $f[0];
        my $Cis_or_Trans= $f[13];
        $Cis_or_Trans =~s/_.*+//g;
        my $Effct_sise = $f[10];
        my $Zscore= "NA";
        unless ($Pos =~/e/){
            my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue";
            print $O1 "$output\n";

            #-----------------------构成vcf文件，用于转hg19;
            my $variant_info = "$Chr\t$Pos\t$Variant_id\t$Ref\t$Alt\t.\t.\t.";
            unless (exists $hash3{$variant_info}){
                $hash3{$variant_info} =1;
                print $O2 "$variant_info\n";
            }
        }
    }
}

while(<$I2>){
    chomp;
    my @f =split/\,/;
    unless(/^SNP/){
        my $Variant_id = $f[0];
        my $Rs_id = "NA";
        my $Chr=$f[1];
        my $Pos= $f[2];
        $Pos =~ s/\.//g;
        my $Ref= $f[7];
        my $Alt= $f[8];
        # print "$Ref\t$Alt\n";
        my $Maf= "NA";
        my $Gene_Version= "Hg18";
        my $Gene_id= "NA";
        my $Slope="NA";
        my $Slope_se= $f[11];
        my $P_value="NA";
        my $Q_value= "NA";
        my $Beta="NA";
        my $Tissue="lymphoblastoid cell lines";
        my $Gene_name= $f[5];
        my $Cis_or_Trans= $f[13];
        $Cis_or_Trans =~s/_.*+//g;
        # print "$Cis_or_Trans\n";
        my $Effct_sise = $f[10];
        my $Zscore= "NA";
        my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effct_sise\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue";
        print $O1 "$output\n";

        #-----------------------构成vcf文件，用于转hg19;
        my $variant_info = "$Chr\t$Pos\t$Variant_id\t$Ref\t$Alt\t.\t.\t.";
        unless (exists $hash3{$variant_info}){
            $hash3{$variant_info} =1;
            print $O2 "$variant_info\n";
        }
    }
}


close $O1;
