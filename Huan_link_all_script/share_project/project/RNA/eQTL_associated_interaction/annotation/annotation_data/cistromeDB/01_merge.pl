#将./Human_FACTOR/human_factor, ./HISTONE_MARK_AND_VARIANT/human_hm, ./Human_CHROMATIN_Accessibility/human_ca下的文件merge到一起"${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz"

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my @dirs = ("./Human_FACTOR/human_factor","./HISTONE_MARK_AND_VARIANT/human_hm","./Human_CHROMATIN_Accessibility/human_ca");
# my @dirs = ("./Human_CHROMATIN_Accessibility/human_ca");
foreach my $dir(@dirs){
    print "$dir\n";
    opendir (DIR, $dir) or die "can't open the directory!";
    my @files = readdir DIR;
    my $output_dir = $dir;
    $output_dir =~s/human.*+//g;
    # print "$output_dir\n";
    my $fo1 = "${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "${output_dir}/merge_pos_info_narrow_peak.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    print $O1 "chr\tstart\tend\tname\tfile_name\n";
    # print $O2 "chr\tstart\tend\n";
    my %hash1;
    my %hash2;
    foreach my $file(@files){
        # print "$file\n";
        my $new_file ="$dir/$file";
        my $f1 = $new_file;
        # # print "$f1\n";
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        print "$f1\n";
        while(<$I1>)
        {
            chomp;
            my @f = split/\s+/;
            my $chr = $f[0];
            my $start =$f[1];
            my $end =$f[2];
            my $name = $f[3];
            my $output1 = "$chr\t$start\t$end\t$name\t$file";
            my $output2 = "$chr\t$start\t$end";
            unless(exists $hash1{$output1}){
                $hash1{$output1}=1;
                print $O1 "$output1\n";
            }
            unless(exists $hash2{$output2}){
                $hash2{$output2}=1;
                print $O2 "$output2\n";
            }
        }
    }
}