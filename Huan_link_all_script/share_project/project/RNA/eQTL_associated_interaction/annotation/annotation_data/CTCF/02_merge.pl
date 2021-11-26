#将data/下面的.gz文件合在一起，并去重得02_all_unique_merge_ctcf_binding_site_narrow_peak.bed.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $dir = "./data";

print "$dir\n";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR;

# print "$output_dir\n";
my $fo1 = "02_all_unique_merge_ctcf_binding_site_narrow_peak.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "chr\tstart\tend\n";
my %hash1;
my %hash2;
foreach my $file(@files){
    # print "$file\n";
    unless($file =~/download/){
        my $new_file ="$dir/$file";
        my $f1 = $new_file;
        print "$f1\n";
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        while(<$I1>)
        {
            chomp;
            unless(/^chr/){
                my @f = split/\s+/;
                my $chr = $f[0];
                my $start =$f[1];
                my $end =$f[2];
                my $output1 = "$chr\t$start\t$end";
                unless(exists $hash1{$output1}){
                    $hash1{$output1}=1;
                    print $O1 "$output1\n";
                }
            }
        }
    }
}