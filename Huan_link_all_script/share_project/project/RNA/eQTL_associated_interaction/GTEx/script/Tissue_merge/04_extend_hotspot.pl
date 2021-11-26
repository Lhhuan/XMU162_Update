#用"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz"对"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz"左边扩9个SNP，右边扩8个SNP得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend.bed.gz",对其进行排序得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;
my $tissue= "Tissue_merge";

my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo3 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;


# print "$tissue\tstart\n";
my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda = $f[0];
        my $pos = $f[1];
        my $chr= $f[2];
        $chr = "chr${chr}";
        $hash1{"$chr\t$pos"}=$.;
        $hash2{$.}="$chr\t$pos";
        # print "$.\t$_\n";
    }
}

print "2 start\n";
while(<$I2>)
{
    chomp;
    my @f=split/\t/;
    my $chr = $f[0];
    my $start = $f[1];
    my $end = $f[2];
    my $true_end= $end -1;
    my $start_k= "$chr\t$start";
    my $end_k= "$chr\t$true_end";
    if (exists $hash1{$start_k}){
        #-----------------------------------new start
        my $start_line = $hash1{$start_k};
        my $new_start_line = $start_line-9; #左扩9，右扩8;
        my $new_start_line_v = $hash2{$new_start_line}; 
        my @t1 =split/\t/,$new_start_line_v;
        my $ns_chr = $t1[0];
        my $ns_pos = $t1[1];
        #------------------------------------new end
        my $end_line = $hash1{$end_k};
        my $new_end_line = $end_line +8;
        my $new_end_line_v = $hash2{$new_end_line}; 
        my @t2 = split/\t/,$new_end_line_v;
        my $ne_chr=$t2[0];
        my $ne_pos = $t2[1];
        if($chr eq $ns_chr && $chr eq $ne_chr){ #这里不判断也可以，因为chr的前9个和后8个pos不用于nhp计算lambda
            my $new_start = $ns_pos;
            my $new_end = $ne_pos+1;
            print $O3 "$chr\t$new_start\t$new_end\n";
        }
    }
}


close($O3);
system "zless $fo3 |sort -k1,1 -k2,2n |gzip >/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted.bed.gz"