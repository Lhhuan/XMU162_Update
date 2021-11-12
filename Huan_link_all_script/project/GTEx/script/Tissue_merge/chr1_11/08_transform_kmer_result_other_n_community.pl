
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


# my $i

# my @aa = (6..10);
my @aa = (3..9);
# my @aa = (3,4);
foreach my $i(@aa){
    my $f1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/communities_${i}.csv";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/communities_${i}.bed.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

    my %hash1;

    while(<$I1>)
    {
        chomp;
        if(/^>/){
            my @f = split/,/;
            my $p = $f[0];
            my $class =$f[1];
            $p=~s/>//g;
            # print "$p\n";
            my @ps= split/:/,$p;
            my $chr = $ps[0];
            my @ses = split/-/,$ps[1];
            my $start = $ses[0];
            my $end = $ses[1];
            print $O1 "$chr\t$start\t$end\t$class\n";

        }
    }

    close($I1);
    close($O1);


    my @factors = ("TFBS","CHROMATIN_Accessibility","CTCF","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3");
    my $file_name = "Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz";
    my $factor_dir = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/annotation";

    my $output_dir ="/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/factor_anno/${i}_community";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    foreach my $factor(@factors){
        my $factor_file = "$factor_dir/${factor}_${file_name}";
        my $command = "bedtools intersect -a $fo1  -b $factor_file -wa |sort -u |gzip >${output_dir}/${factor}.bed.gz" ;
        print "$factor_file\n";
        system $command;
    }
}