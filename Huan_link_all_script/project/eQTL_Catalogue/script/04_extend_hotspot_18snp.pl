#用"../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz"对"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz"扩展，只对<18 个snp的hotspot进行扩展，以hotspot中心开始，扩到18个snp,得"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18snp.bed.gz",对其进行排序得"../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted.bed.gz"
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

my $f1 = "../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "./count_snp/output/02_count_snp_in_hotspot_sorted.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo3 = "../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18snp.bed.gz";
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
        my $pos_bed  =$pos-1;
        my $chr= $f[2];
        $chr = "chr${chr}";
        $hash1{"$chr\t$pos_bed"}=$.;
        $hash2{$.}="$chr\t$pos_bed";
        # print "$.\t$_\n";
    }
}

print "2 start\n";
while(<$I2>)
{
    chomp;
    unless(/SNP_number/){
        my @f=split/\t/;
        my $chr = $f[0];
        my $start = $f[1];
        my $end = $f[2];
        my $SNP_number =$f[3];
        my $center_snp=$f[4];
        if($SNP_number >=18){
            print $O3 "$_\n";
        }
        else{
            my $k = "$chr\t$center_snp";
            if(exists $hash1{$k}){
                my $center_line=$hash1{$k};
                #------------------start
                my $new_start_line=$center_line-9;
                my $new_start_line_v = $hash2{$new_start_line}; 
                my @t1 =split/\t/,$new_start_line_v;
                my $ns_chr = $t1[0];
                my $ns_pos = $t1[1];
                #-------------------end
                my $new_end_line = $center_line+8;
                my $new_end_line_v = $hash2{$new_end_line}; 
                my @t2 = split/\t/,$new_end_line_v;
                my $ne_chr=$t2[0];
                my $ne_pos = $t2[1];
                my $new_start = $ns_pos-1+1; #防止科学计数法 #因为前面将输入改成0-based
                my $new_end = $ne_pos+1;#防止科学计数法
                print $O3 "$chr\t$new_start\t$new_end\t$SNP_number\t$center_snp\n";
            }
        }
    }   
}


close($O3);
system "zless $fo3 |sort -k1,1 -k2,2n |gzip > ../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted.bed.gz"