 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 ="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/cis_eQTLs_all_re.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

#---------------
my $fo1 = "../output/021_LAML_cis_eQTL_sig_egene.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "SNP_chr\tSNP_pos\tPvalue\tegene\n";

my $tissue = "LAML";
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
        my $gene = $f[5];
        my $Pvalue =$f[-1];
        if($cancer_type eq $tissue){
            print "$cancer_type\n";
            if($Pvalue <5e-8){
                # my $k = "$SNP_chr\t$SNP_pos";
                  print $O1 "$SNP_chr\t$SNP_pos\t$Pvalue\t$gene\n";
            }
        }
    }
}
