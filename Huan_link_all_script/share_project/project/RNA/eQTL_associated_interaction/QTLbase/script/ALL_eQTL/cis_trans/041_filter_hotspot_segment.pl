#将../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot.bed.gz 长度>=6的片段，得
#../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot_length_more_than6.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @QTLs =("eQTL");
my $pm = Parallel::ForkManager->new(7); ## 设置最大的线程数目
foreach my $QTL(@QTLs){
    # my @interval = (18);
    # my @interval = (15);
    my @interval = (6,7,8,9,12,15,18);
    my @types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
    foreach my $type(@types){
        foreach my $j(@interval){
            my $pid = $pm->start and next; #开始多线程
            print "$j\t$type\n";
            my $f1 = "../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot.bed.gz";
            # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
            open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
            my $fo1 = "../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot_length_more_than6.bed.gz";
            open my $O1, "| gzip >$fo1" or die $!;
            my %hash1;
            while(<$I1>)
            {
                chomp;
                my @f = split/\t/;
                my $chr = $f[0];
                my $start = $f[1];
                my $end= $f[2];
                my $length = $end -$start;
                if ($length >=6){
                    print $O1 "$_\n";
                }
            }
            # ---------------------------
            

            $pm->finish;  #多线程结束  
        }
    }
}
