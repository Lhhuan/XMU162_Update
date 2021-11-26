#对"E098-H3K27ac.narrowPeak.gz","E098-H3K4me1.narrowPeak.gz" 进行normal和排序得${fo1}sorted.gz
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
# use Parallel::ForkManager;
my @files = ("E098-H3K27ac.narrowPeak.gz","E098-H3K4me1.narrowPeak.gz");
foreach my $file(@files){
    my $f1 = $file;
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    
    my $out= $file;
    $out =~s/\.gz//g;
    my $fo1 = "01_normal_${out}";
    # open my $O1, "| gzip >$fo1" or die $!;
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    my %hash1;
    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $chr= $f[0];
        # $chr =~s/chr//g;
        my $start= $f[1];
        my $end= $f[2];
        unless($chr=~/chrY|chrX/){
            print $O1 "$chr\t$start\t$end\n";
        }
    }
    my $command1 = "cat $fo1 | sort -k1,1 -k2,2n  | gzip > ${fo1}sorted.gz";
    print "$command1\n";
    system $command1;
    
}
