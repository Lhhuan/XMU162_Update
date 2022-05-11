#利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz"为 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290_sorted.bed(../output/01_hotspot_length_small_than_2290_sorted_egene.bed.gz)寻找egene,得 "../output/01_hotspot_length_small_than_2290_sorted_egene_h.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;



my $f1 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
my $fo2 = "../output/01_hotspot_egene_h.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "Chr\tstart\tend\tegene\thotspot\n";


while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $chr=$f[0];
    if($chr=~/^chr1$/ | $chr=~/^chr22$/){
        print "$chr\n";
        my $hotspot = join("_",@f[0..2]);
        print $O2 "$_\t$hotspot\n";
    }

}
