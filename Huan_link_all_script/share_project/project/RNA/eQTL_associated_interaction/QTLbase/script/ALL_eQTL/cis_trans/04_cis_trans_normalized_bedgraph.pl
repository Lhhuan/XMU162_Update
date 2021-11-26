#将../output/ALL_eQTL/NHPoisson_emplambda_interval_${i}cutoff_7.3_all_eQTL.txt.gz 中的 chr1中的数据提出来得../output/ALL_eQTL/chr1/NHPoisson_emplambda_interval_${i}cutoff_7.3_all_eQTL.bedgraph
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my @types = ("cis_1MB","cis_10MB","trans_10MB","trans_1MB");
# my @number = (10..18);
my @number = (10..20);
my $pm = Parallel::ForkManager->new(40); ## 设置最大的线程数目
foreach my $type(@types){
    foreach my $i(@number){
        my $pid = $pm->start and next; #开始多线程
        my $f1 = "../../../output/ALL_eQTL/cis_trans/NHP/${type}/NHPoisson_emplambda_interval_${i}_cutoff_7.3_${type}_eQTL.txt.gz";
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

        # my $fo1 = "../output/ALL_eQTL/chr1/NHPoisson_emplambda_interval_${i}cutoff_7.3_all_eQTL.bedgraph";
        my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/chr1_2/${type}/NHPoisson_emplambda_interval_${i}cutoff_7.3_${type}_eQTL.bedgraph";
        open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        # open my $O1, "| gzip >$fo1" or die $!;


        while(<$I1>)
        {
            chomp;
            unless(/^emplambda|NA/){
                my @f = split/\t/;
                my $emplambda=$f[0];
                my $POS =$f[1];
                my $chr =$f[2]; 
                if ($chr eq 1 || $chr eq 2){
                    # print "$chr\n";
                    my $start = $POS;
                    my $end = $start +1;
                    $start =$end -1; 
                    # print "$i\n";
                    # my $CHR = "chr${chr}";
                    print $O1 "$chr\t$start\t$end\t$emplambda\n";
                }
            }
        }
        $pm->finish;  #多线程结束
    }
}