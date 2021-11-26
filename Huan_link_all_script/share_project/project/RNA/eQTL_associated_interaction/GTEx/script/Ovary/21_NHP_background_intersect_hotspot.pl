#将"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz” 转换为bed得 "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz,排序得$sort_fo1，用$sort_fo1 和"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz bedtools intersect得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;


my $tissue= "Ovary";
my %hash1;
my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;

while(<$I1>)
{
    chomp;
    unless(/^SNP_chr/){
        my @f = split/\t/;
        my $SNP_chr = $f[0];
        my $SNP_pos = $f[1];
        my $Pvalue =$f[2];
        my $end = $SNP_pos +1;
        my $chr = "chr${SNP_chr}";
        print $O1 "$chr\t$SNP_pos\t$end\t$Pvalue\n";
        
    }
}

close($I1);
close($O1);
my $sort_fo1 = $fo1;
$sort_fo1 =~ s/_Completion/_Completion_sorted/g;
my $command1 = "zless $fo1 | sort -k1,1 -k2,2n |gzip > $sort_fo1";
system $command1;
my @cutoffs;
for (my $i=0.05;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
}
my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
foreach my $cutoff(@cutoffs){
    my $pid = $pm->start and next; #开始多线程
    my $hotspot = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $output_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    my $out_file = "${output_dir}/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz";
    my $command2 = "bedtools intersect -a $sort_fo1 -b $hotspot -wa -wb | gzip > $out_file";
    system $command2;
    print "$cutoff\n";
    $pm->finish;  #多线程结束
}

