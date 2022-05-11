#将./enhancer_gene 下的文件合并起来得 02_merge_enhancer_target_sample.bed.gz, sorted  02_merge_enhancer_target_sample_sorted.bed.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
my $dir = "./enhancer_gene";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR;
#---------------------------------output option
my $fo1 = "02_merge_enhancer_target_sample.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
#--------对dir 下all file 进行判断
foreach my $file(@files){
    if ( $file =~ /txt/) {
        my $cell_line= $file;
        $cell_line =~s/\_EP.*//g;
        print "$cell_line\n";
        my $f3 = "$dir/$file";
        open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
        # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
        # print "$f3\n";
        while(<$I3>)
        {
            chomp;
            my @q=split/\s+/;
            my $pred_score = $q[1];
            my @f = split/\_/;
            my $pos= $f[0];
            my $gene =$f[1];
            my @t=split/:/,$pos;
            my $chr = $t[0];
            my @ff =split/-/,$t[1];
            my $start =$ff[0];
            my $end =$ff[1];

            $gene =~ s/\$.*//g;
            print $O1 "$chr\t$start\t$end\t$gene\t$pred_score\t$cell_line\n";
        }
    }
}
close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > 02_merge_enhancer_target_sample_sorted.bed.gz";
