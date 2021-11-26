#统计每个community中各个部分显著motif"${output_dir}/knowResult_merge.txt"， overlap情况得"${output_dir}/knowResult_merge_class_merge.txt",overlap ratio 文件"./6kmer/community_motif_overlap_ratio.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;



# my $fo3 = "community_motif_overlap_ratio.txt";
# open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
# print $O3 "Number_of_all_motif\tNumber_of_overlap_motif\toverlap_ratio\tcommunity\n";

my %hash1;
my $f1 ="/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/knowResult_merge.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/knowResult_specific.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Motif_Name\tConsensus\tcommunity\n";
while(<$I1>)
{
    chomp;
    unless(/^Motif/){
        my @f= split/\t/;    
        my $motif =$f[0] ;
        my $Consensus =$f[1];
        my $community = $f[3];
        my $k = "$motif\t$Consensus";
        push @{$hash1{$k}},$community;
    }

} 


foreach my $k(sort keys %hash1){
    my @vs = @{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs;
    my $number =@vs;
    if ($number <2){
        print $O1 "$k\t$vs[0]\n";
        print "@vs\n";
    }
}    
