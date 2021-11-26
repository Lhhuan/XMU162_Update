
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Whole_Blood";
my $state =15;
my $cutoff2 = "0";

my $dir = "/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2";
my $pm = Parallel::ForkManager->new(40); ## 设置最大的线程数目
for (my $i=1;$i<10001;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $name = "${state}_state_resemble_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    my $f1 = "$dir/$name";

    my $command = "rm $f1";
    system $command;
    print "$i\n";
    $pm->finish;  #多线程结束
}





