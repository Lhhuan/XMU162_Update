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
my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
for (my $i=1;$i<1001;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $input_file = "../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/original_random/${cutoff}/sorted_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    my $input_file_base_name = basename($input_file);
    my $output_dir = "../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/original_random/${cutoff}_filter_6"; 
    # mkdir $PMID;
    # #------------
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    #------------

    my $f1 =$input_file;
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    #----------------------
    my $fo1 = "${output_dir}/$input_file_base_name";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    while(<$I1>)
    {
        chomp;
        my @f= split/\t/;
        my $SNP_chr =$f[0];
        my $start =$f[1];
        my $end =$f[2];
        my $length = $end-$start;
        if($length >5){
            print $O1 "$_\n";
        }
    }
    print "$i\t$group\t$cutoff\n";
    $pm->finish;  #多线程结束
}
      