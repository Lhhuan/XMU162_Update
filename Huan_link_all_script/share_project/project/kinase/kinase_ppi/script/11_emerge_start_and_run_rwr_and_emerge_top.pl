 ###../output/08_all_gene_from_jintao_in_the_network.txt 每个基因symbol作为成RWR的起点，并基于../output/network_symbol_used_rwr.txt
#run rwr，并把结果存在../output/RWR_a_g_JT/rwr_result/文件夹下面
#取 rwr result的top 存在../output/RWR_a_g_JT/rwr_result_top/文件夹下面
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use File::Copy;



my $f1 ="../output/08_all_gene_from_jintao_in_the_network.txt";#输入的是start
 open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
 my %hash1;

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^symbol/){
        # my $query = $f[0];
        my $symbol = $f[0];
        my $fo1 ="../output/RWR_a_g_JT/start/${symbol}.txt"; 
        open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        print $O1 "$symbol\n";
        close $O1 or warn "$01 : failed to close output file '$fo1' : $!\n";
        print STDERR "$symbol\n";
        system "python run_walker.py ../output/network_symbol_used_rwr.txt ../output/RWR_a_g_JT/start/${symbol}.txt > ../output/RWR_a_g_JT/rwr_result/${symbol}.txt";
        system "cat ../output/RWR_a_g_JT/rwr_result/${symbol}.txt | sort -k2,2rg >../output/RWR_a_g_JT/rwr_result/${symbol}_sorted.txt";
        # 13852
        my $line = 0.01*13852;
        my $line2  = sprintf "%.f", $line; # 这个是四舍五入取整
        system "head -n $line2 ../output/RWR_a_g_JT/rwr_result/${symbol}_sorted.txt > ../output/RWR_a_g_JT/rwr_result_top/${symbol}.txt"; #取top
    }
}

