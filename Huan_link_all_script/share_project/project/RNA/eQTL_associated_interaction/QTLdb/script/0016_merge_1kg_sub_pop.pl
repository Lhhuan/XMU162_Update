#将/share/data0/1kg_phase3_v5_hg19/五个人种的all_SNP merge在一起得../output/all_1kg_phase3_v5_hg19_snp.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "../output/all_1kg_phase3_v5_hg19_snp.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "CHROM\tPOS\tID\tref\talt\n";
#------------------------------------
my @files = ("/share/data0/1kg_phase3_v5_hg19/AFR/1kg.phase3.v5.shapeit2.afr.hg19.chr9.vcf.gz","/share/data0/1kg_phase3_v5_hg19/SAS/1kg.phase3.v5.shapeit2.sas.hg19.all.SNPs.vcf.gz","/share/data0/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz","/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz",
"/share/data0/1kg_phase3_v5_hg19/AMR/1kg.phase3.v5.shapeit2.amr.hg19.all.SNPs.vcf.gz");
my %hash1;
foreach my $file(@files){
    # print "$file\n";
    my @t= split/\//,$file;
    my $pop = $t[4];
    # print "$pop\n";
    my $f1 = $file;
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    

    while(<$I1>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\s+/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $ID =$f[2]; #snp chr
            my $ref =$f[3];
            my $alt = $f[4];
            my $output = "$CHROM\t$POS\t$ID\t$ref\t$alt";
            unless (exists $hash1{$output}){
                $hash1{$output}=1;
                print $O1 "$output\n";
            }
        }
    }
}