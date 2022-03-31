#将"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz" split 成chr1-22,得"../output/chr_split/01_merge_all_tissue_cis_eQTL_1kg_Completion_chri.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

# my $pm = Parallel::ForkManager->new(30); ## 设置最大的线程数目
# foreach my $i(1..22){
    # my $pid = $pm->start and next; #开始多线程
    my $i=10;
    print "$i\tstart\n";
    my $fo2 = "../output/chr_split/01_merge_all_tissue_cis_eQTL_1kg_Completion_chr${i}.txt.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    print $O2 "SNP_chr\tSNP_pos\tPvalue\n";

    my $f1 = "../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        unless(/^SNP_chr/){
            my $chr =$f[0];
            if($chr == $i){
                print $O2 "$_\n";
            }
        }
    }
    close($I1);
    close($O2);
    print "$i\tend\n";
    # $pm->finish;  #多线程结束

# }
