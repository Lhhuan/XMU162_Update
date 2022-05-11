#将"../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz" split 成chr1-22,得"../output/02_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted_chr${i}.bed.gz"
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $i=11;
print "$i\tstart\n";
my $fo2 = "../../output/02_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted_chr${i}.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# print $O2 "SNP_chr\tSNP_pos\tPvalue\n";

my $f1 = "../../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    unless(/^SNP_chr/){
        my $chr =$f[0];
        my $q="chr${i}";
        if($chr eq $q){
            print $O2 "$_\n";
        }
    }
}
close($I1);
close($O2);
print "$i\tend\n";

