 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $tissue = "LAML";

print "$tissue\tstart\n";
my $f1 ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/cis_eQTLs_all_re.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $f2 = "/share/data0/QTLbase/huan/1kg_phase3_v5_hg19/1kg.phase3.v5.shapeit2.hg19.all.vcf.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件


#---------------
my $fo1 = "../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "SNP_chr\tSNP_pos\tPvalue\n";


my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^cancer_type/){
        my @f = split/\t/;
        my $cancer_type =$f[0];
        my $SNP_chr =$f[2];
        $SNP_chr =~ s/chr//g;
        my $SNP_pos =$f[3];
        my $Pvalue =$f[-1];
        if($cancer_type eq $tissue){
            # print "$cancer_type\n";
            my $k = "$SNP_chr\t$SNP_pos";
            $hash1{$k}=1;
            print $O1 "$k\t$Pvalue\n";
        # print "$k\t$Pvalue\n";
        }
    }
}


while(<$I2>)
{
    chomp;
    unless(/^#/){
        my @f = split/\t/;
        my $CHROM =$f[0];
        my $POS =$f[1]; 
        my $pvalue = 0.05;
        my $k = "$CHROM\t$POS";
        unless (exists $hash1{$k}){
            print $O1 "$k\t$pvalue\n";
            $hash1{$k}=1;
        }
    }
}
