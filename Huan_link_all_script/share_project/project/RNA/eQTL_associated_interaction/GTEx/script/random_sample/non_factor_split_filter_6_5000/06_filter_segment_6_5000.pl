#interval_18 时，对../../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
#extend 100bp
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my @intervals = (18,15,12,9,8,7,6);
my @cutoffs = ();
for (my $i=0.1;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
my @groups = ("hotspot");
my @types = ("non_factor_split","factor");
my @tissues = ("Whole_Blood","Lung");
my @numbers= (1..10);

my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
foreach my $tissue(@tissues){
    foreach my $type(@types){
        foreach my $group(@groups){
            foreach my $number(@numbers){
                foreach my $cutoff(@cutoffs){
                    my $pid = $pm->start and next; #开始多线程
                    my $input_file = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/${group}/${number}/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
                    # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
                    my $input_file_base_name = basename($input_file);
                    my $output_dir = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/${group}/${number}/interval_18_filter/6_5000"; 
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
                        if($length >5 && $length <5001){
                            print $O1 "$_\n";
                        }
                    }
                    print "$type\t$group\t$cutoff\n";
                    $pm->finish;  #多线程结束
                }
            } 
        }
    }
}

