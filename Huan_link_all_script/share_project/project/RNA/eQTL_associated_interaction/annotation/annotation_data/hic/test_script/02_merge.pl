#将hg19/下面的.loops文件合在一起，由这些文件判断，该文件是1-based，改为0-based,判断cis-trans_interaction 得02_all_hic_loops_tran_cis_1.bed.gz (interaction1_interaction2),
#得02_all_hic_loops_tran_cis_2.bed.gz， 得不去重的文件得02_all_hic_loops_source.bed.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $dir = "./hg19";

print "$dir\n";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR;

my $fo1 = "02_all_hic_loops_source.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "02_all_hic_loops_tran_cis_1.bed.gz";#原来的顺序interaction1_interaction2
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O2, "| gzip >$fo2" or die $!;
my $fo3 = "02_all_hic_loops_tran_cis_2.bed.gz"; #把interaction2_interaction1
open my $O3, "| gzip >$fo3" or die $!;
print $O1 "chr1\tstart1\tend1\tchr2\tstart2\tend2\tsource_file_name\n";
print $O2 "chr1\tstart1\tend1\tchr2\tstart2\tend2\t\tcis_trans_1MB\tcis_trans_10MB\n";
print $O3 "chr2\tstart2\tend2\tchr1\tstart1\tend1\t\tcis_trans_1MB\tcis_trans_10MB\n";
my (%hash1,%hash2,%hash3,%hash4);
foreach my $file(@files){
    # print "$file\n";
    my $new_file ="$dir/$file";
    my $f1 = $new_file;
    print "$f1\n";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    while(<$I1>)
    {
        chomp;
        my @f = split/\s+/;
        my $chr1 = $f[0];
        my $start1 =$f[1];
        my $end1 =$f[2];
        my $chr2 = $f[3];
        my $start2 =$f[4];
        my $end2 =$f[5];
        my $new_start1 = $start1 -1; #转化为0-based
        my $new_start2 = $start2 -1; #转化为0-based
        my $output2 = "$chr1\t$new_start1\t$end1\t$chr2\t$new_start2\t$end2";
        my $output = "$output2\t$file";
        print $O1 "$output\n";
        my $output3 = "$chr2\t$new_start2\t$end2\t$chr1\t$new_start1\t$end1";
        my $different = $new_start1 -$new_start2;
        print $O2 "$output2\t";
        print $O3 "$output3\t";
        if($chr1 eq $chr2 && $different <=1000000){ #1MB_cis
            print $O2 "cis\t";
            print $O3 "cis\t";
        }
        else{ ##1MB_trans
            print $O2 "trans\t";
            print $O3 "trans\t";
        }
        if ($chr1 eq $chr2 && $different <=10000000){ #10MB_cis
            print $O2 "cis\n";
            print $O3 "cis\n";
        }
        else{##10MB_trans
            print $O2 "trans\n";
            print $O3 "trans\n";
        }
    }
}