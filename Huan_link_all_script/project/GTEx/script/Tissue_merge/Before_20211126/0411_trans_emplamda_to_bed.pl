###../../output/${tissue}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${tissue}.txt.gz 时得不同cutoff下的hotspot(segment),"../../output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;


my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo3 = "/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda = $f[0];
        my $pos = $f[1];
        my $chr= $f[2];
        my $start = $pos;
        my $end =$start +1;
        $start = $end-1;
        unless($emplambda =~/NA/){
           print $O3 "chr${chr}\t$start\t$end\t$emplambda\n";
        }
    }
}
close($O3);
system "zless $fo3 |sort -k1,1 -k2,2n |gzip >/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge_sorted.bed.gz"