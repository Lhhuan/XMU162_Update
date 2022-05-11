#将./enhancer 下的文件合并起来得 01_merge_enahcner_sample.bed.gz 

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $dir = "./enhancer";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR;
#---------------------------------output option
my $fo1 = "01_merge_enahcner_sample.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
#--------对dir 下all file 进行判断
foreach my $file(@files){
    if ( $file =~ /bed/) {
        my $cell_line= $file;
        $cell_line =~s/\.bed//g;
        print "$cell_line\n";
        my $f3 = "$dir/$file";
        open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
        # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
        # print "$f3\n";
        while(<$I3>)
        {
            chomp;
            my @f = split/\s+/;
            my $chr = $f[0];
            my $start =$f[1];
            my $end =$f[2];
            print $O1 "$chr\t$start\t$end\t$cell_line\n";
        }
    }
}
close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > 01_merge_enahcner_sample_sorted.bed.gz";
