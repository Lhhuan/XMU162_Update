#将../raw_data/2012-12-21-CisAssociationsProbeLevelFDR0.5.txt 和../raw_data/2012-12-21-TransEQTLsFDR0.5.txt 进行初步normal, 得gene name 和部分RSid从对应的v8.egenes.txt.gz中获得，得文件../output/01_normal_format.txt.gz
#并得unique的varint POS的vcf文件../output/01_unique_variant-id_position_hg18.vcf.gz
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
my @files = ("2012-12-21-CisAssociationsProbeLevelFDR0.5.txt","2012-12-21-TransEQTLsFDR0.5.txt");
my(%hash1,%hash2,%hash3,%hash4);
foreach my $file(@files){
    my $f1 = "../raw_data/${file}";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
    # print "$f1\n";
    while(<$I1>){
        chomp;
        my @f =split/\t/;
        unless(/^PValue/){
            my $Variant_id = $f[1];
            my $Rs_id = $f[1];
            my $Chr=$f[2];
            my $Pos= $f[3];
            my $SNPType = $f[8];
            my @t =split/\//,$SNPType;
            my $Ref= $t[0];
            my $Alt= $t[1];
            # print "$Ref\t$Alt\n";
            my $Maf= "NA";
            my $Gene_Version= "Hg18";
            my $Gene_id= "NA";
            my $Slope="NA";
            my $Slope_se= "NA";
            my $P_value=$f[0];
            my $Q_value= "NA";
            my $Beta= "NA";
            my $Tissue="Blood" ;
            my $Effect_size = "NA";
            my $Zscore= "NA";
            if ($file =~/Cis/){ #两个文件gene name 的位置不同
                my $Gene_name= $f[-2] ;
                my $Cis_or_Trans= "cis";
                my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effect_size\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue";
                    print $O1 "$output\n";
            }
            else{
                my $Gene_name= $f[16];
                my $Cis_or_Trans= "trans";
                my $output= "$Variant_id\t$Rs_id\t$Chr\t$Pos\t$Ref\t$Alt\t$Maf\t$Gene_Version\t$Gene_id\t$Gene_name\t$Slope\t$Slope_se\t$Zscore\t$Effect_size\t$P_value\t$Q_value\t$Beta\t$Cis_or_Trans\t$Tissue";
                    print $O1 "$output\n";
            }
            #-----------------------构成vcf文件，用于转hg19;
            my $variant_info = "$Chr\t$Pos\t$Variant_id\t$Ref\t$Alt\t.\t.\t.";
            unless (exists $hash3{$variant_info}){
                $hash3{$variant_info} =1;
                print $O2 "$variant_info\n";
            }
        }
    }
    close $I1;
}
close $O1;
