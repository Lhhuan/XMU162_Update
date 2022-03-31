#利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz"为 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed(../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz)寻找egene,得 "../output/01_hotspot_length_small_than_2290_sorted_egene_h.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $fo1 = "../output/01_merge_all_tissue_cis_eQTL_gene.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Chr\tStart\tEnd\tegene\thotspot\ttissue\n";

my $f2 = "/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19/Brain_Spinal_cord_cervical_c-1.v8.signif_variant_gene_pairs.txt.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my (%hash1,%hash2);

while(<$I2>)
{
    chomp;
    my @f =split/\t/;
    unless(/^variant_id/){
        my $SNP_chr = $f[1];
        my $SNP_pos =$f[2]; #0-based
        my $Pvalue=$f[-6];
        my $gene_id=$f[3];
        $gene_id =~ s/\..*+//g;
        my $start = $SNP_pos;
        my $end =$SNP_pos +1;
        my $chr="chr${SNP_chr}";
        # if($Pvalue <5E-8){
            # print "$f[3]\t$gene_id\n";
        my $output= "$chr\t$start\t$end\t$gene_id";
        unless(exists $hash2{$output}){
            $hash2{$output}=1;
            print $O1 "$output\n";
        }
    }
}

close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n| gzip >../output/01_merge_all_tissue_cis_eQTL_gene_sorted.bed.gz";


system "bedtools intersect -a /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed -b ../output/01_merge_all_tissue_cis_eQTL_gene_sorted.bed.gz -wa -wb |cut -f1-3,7|sort -u |gzip >../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz";

my $f1 = "../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
my $fo2 = "../output/01_hotspot_length_small_than_2290_sorted_egene_h.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "Chr\tstart\tend\tegene\thotspot\n";


while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $hotspot = join("_",@f[0..2]);
    print $O2 "$_\t$hotspot\n";
}
