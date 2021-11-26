#将$input_file 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group ="hotspot";
my $tissue = "Whole_Blood";
my $cutoff2 = "original_genome";
chdir "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/$cutoff2/fa";
$ENV{'cutoff'}  = $cutoff; #设置环境变量
$ENV{'group'} = $group ;
$ENV{'tissue'} = $tissue ;
my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
for (my $i=1;$i<1001;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $input = "${i}_resemble_extend_${tissue}_segment_${group}_cutoff_${cutoff}.fa";
    $ENV{'i'} = $i;
    $ENV{'input'} = $input ;
    my $command1 = "bash /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/get_kmer.sh";
    system $command1;
    print "$i\n";
    $pm->finish;  #多线程结束
}
      