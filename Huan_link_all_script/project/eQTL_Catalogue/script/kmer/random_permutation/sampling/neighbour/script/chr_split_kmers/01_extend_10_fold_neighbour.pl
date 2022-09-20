#利用"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data/need_download_tabix_ftp_paths.tsv" 部分文件，提取pos，p,gene得"../output/01_merge_all_tissue_cis_eQTL_eur_egene.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $fo1 = "../output/01_extend_10_fold_neighbour.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $f1 = "/share/Projects/huanhuan/ref_data/UCSC/hg38/hg38.chrom1_22_sizes_sorted.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../../../../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my %hash1;
while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr  =$f[0];
    my $length = $f[1];
    $hash1{$chr}=$length;
}

while(<$I2>)
{
    chomp;
    my @f = split/\t/;
    my $chr  =$f[0];
    my $start =$f[1];
    my $end=$f[2];
    my $length = $end - $start;
    my $new_start =$start - 10*$length;
    my $new_end= $end + 10*$length;
    my $chr_len = $hash1{$chr};
    if ($new_start >=0 && $new_end <=$chr_len ){
        print $O1 "$chr\t$new_start\t$new_end\n";
    }
    elsif($new_start <0 && $new_end >$chr_len ){
        print $O1 "$chr\t0\t$chr_len\n";
    }
    elsif($new_start <0 && $new_end <=$chr_len){
        print $O1 "$chr\t0\t$new_end\n";
    }
    else{ #$new_start >=0 && $new_end >$chr_len
        print $O1 "$chr\t$new_start\t$chr_len\n";
    }
    
}

close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ../output/01_extend_10_fold_neighbour_sort.bed.gz";
system "bedtools merge -i  ../output/01_extend_10_fold_neighbour_sort.bed.gz |gzip > ../output/01_extend_10_fold_neighbour_sort_merge.bed.gz";
