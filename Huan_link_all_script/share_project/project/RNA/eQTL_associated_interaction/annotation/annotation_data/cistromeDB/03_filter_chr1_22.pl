#过滤出${output_dir}/merge_pos_info_narrow_peak_sorted.bed.gz中的chr1_22,得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz
#将${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz左右各扩500bp得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_500bp.bed.gz
#将${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz左右各扩100bp得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;



my @dirs = ("./Human_FACTOR","./HISTONE_MARK_AND_VARIANT","./Human_CHROMATIN_Accessibility");
# my $pm = Parallel::ForkManager->new(4); ## 设置最大的线程数目
foreach my $dir(@dirs){
    # my $pid = $pm->start and next; #开始多线程
    my $output_dir = $dir;
    print "$dir\n";
    my $fo1 = "${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_500bp.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O2, "| gzip >$fo2" or die $!;
    my $fo3 = "${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O3, "| gzip >$fo3" or die $!;



    my $f1 = "${output_dir}/merge_pos_info_narrow_peak_sorted.bed.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my %hash1;
    my %hash2;

    for (my $i=1;$i<23;$i++){
        my $k = "chr${i}";
        $hash1{$k}=1;
    }


    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $chr = $f[0];
        if (exists $hash1{$chr}){
            print $O1 "$_\n";
            my $chr = $f[0];
            my $start = $f[1];
            my $end = $f[2];
            my $extend_500_start = $start-500;
            my $extend_500_end = $end+500;
            my $extend_100_start = $start-100;
            my $extend_100_end = $end+100; 
            if($extend_500_start >=0){
                print $O2 "$chr\t$extend_500_start\t$extend_500_end\t$start\t$end\n";
            }
            else{
                print $O2 "$chr\t$start\t$end\t$start\t$end\n";
            }
            if($extend_100_start >=0){
                print $O3 "$chr\t$extend_100_start\t$extend_100_end\t$start\t$end\n"; 
            }
            else{
                print $O3 "$chr\t$start\t$end\t$start\t$end\n"; 
            }
        }
    }
    # $pm->finish;  #多线程结束

}