#将../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz 中的QTLbase中../output/merge_QTL_all_QTLtype_pop.txt.gz的数据过滤出来，
#得在qtlbase文件../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz
#得不在qtlbase文件../output/ALL_${xQTL}/1kg_out_QTLbase_NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;

my (%hash1,%hash2,%hash4);
my $f1 = "../output/merge_QTL_all_QTLtype_pop.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件


while(<$I1>)
{
    chomp;
    unless(/SNP_chr/){
        my @f = split/\t/;
        my $snp_chr = $f[0];
        my $snp_pos = $f[1];
        my $QTL_type = $f[8];
        my $k1 = "$snp_chr\t$snp_pos\t$QTL_type";
        my $k2 = "$snp_chr\t$snp_pos";
        $hash1{$k1}=1;
        $hash2{$k2}=1;

        
    }
}


# my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL","cerQTL","lncRNAQTL","metaQTL","miQTL","riboQTL");
my @QTLs = ("eQTL");
# my $xQTL = "eQTL" ;
my $pm = Parallel::ForkManager->new(15); ## 设置最大的线程数目
foreach my $xQTL(@QTLs){
    my $pid = $pm->start and next; #开始多线程
    my $f2 = "../output/ALL_${xQTL}/NHPoisson_emplambda_interval_1000cutoff_7.3_all_${xQTL}.txt.gz";
    # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    my $fo1 = "../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "../output/ALL_${xQTL}/1kg_out_QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_${xQTL}.txt.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    print "start\t$xQTL\n";
    my %hash3;
    while(<$I2>)
    {
        chomp;
        if(/emplambda/){
            print $O1 "$_\n";
            print $O2 "$_\n";
        }
        else{
            my @f = split/\t/;
            my $emplambda =$f[0];
            my $t =$f[1]; 
            my $pos = $t; #snp pos
            my $chr =$f[2]; #snp chr
            my $k1 = "$chr\t$pos\t$xQTL";
            my $k2 = "$chr\t$pos";
            if ($xQTL =~/\bQTL\b/){ #---all QTL
                if(exists $hash2{$k2}){
                    print $O1 "$_\n";
                }
                else{
                    print $O2 "$_\n";
                }
            }
            else{ #------other QTL
                if (exists $hash1{$k1}){
                    $hash3{$k1}=1;
                    print $O1 "$_\n";
                }
                else{
                    print $O2 "$_\n";
                }
            }
            
        }
    }
    close($I2);
    close($O1);
    close($O2);
    $pm->finish;  #多线程结束
}