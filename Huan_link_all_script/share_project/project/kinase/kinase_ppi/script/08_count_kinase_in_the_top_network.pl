#根据../output/07_final_all_gene_from_jintao_symbol.txt和../output/RWR/rwr_result_top/ 计算kinase in the top, 
#and kinase out of top, isn't no-kinase in the top and isn't no-kinase out of the top,得文件../output/08_fisher_need_data.txt
#得不在网络中的非kinase gene文件../output/08_nonkinase_not_in_the_network.txt 
#得在网络中的非kinase gene文件 ../output/08_nonkinase_in_the_network.txt
#得不在网络中的kinase gene文件../output/08_kinase_not_in_the_network.txt
#得在网络中的kinase gene文件 ../output/08_kinase_in_the_network.txt
#得symbol的$neighbour文件 ../output/08_start_neighbour.txt
#得 在不同start下的，neighbour得event文件 ../output/08_start_neighbour_events.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use File::Copy;

my $f1 ="../output/03_final_all_network_gene_symbol.txt";
my $f2 ="../output/07_final_all_gene_from_jintao_symbol.txt";#
my $f3 ="../output/01_transfrom_kinsae_gene_to_symbol.txt";#输入的是start

open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
my $fo1 ="../output/08_nonkinase_not_in_the_network.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 ="../output/08_nonkinase_in_the_network.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 ="../output/08_kinase_in_the_network.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my $fo4 ="../output/08_kinase_not_in_the_network.txt";
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";
my $fo5 ="../output/08_fisher_need_data.txt";
open my $O5, '>', $fo5 or die "$0 : failed to open output file '$fo5' : $!\n";
my $fo6 ="../output/08_start_neighbour.txt";
open my $O6, '>', $fo6 or die "$0 : failed to open output file '$fo6' : $!\n";
my $fo7 ="../output/08_start_neighbour_events.txt";
open my $O7, '>', $fo7 or die "$0 : failed to open output file '$fo7' : $!\n";
print $O1 "symbol\n";
print $O2 "symbol\n";
print $O3 "symbol\n";
print $O4 "symbol\n";
print $O5 "symbol\tevent_in_the_top_final\tgene_not_event_in_the_top\tevent_out_of_top\tgene_not_event_out_of_the_top\n";
print $O6 "symbol\tneighbour\n";
print $O7 "symbol\tevents_in_neighbour\n";
my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^query/){
        my $query = $f[0];
        my $symbol = $f[1];
        $hash1{$symbol}=1;
    }
}

while(<$I2>)
{
    chomp;
    my @f= split /\t/;
    unless(/^Gene/){
        my $query = $f[0];
        my $symbol = $f[1];
        my $kinase = $f[-1];
        if($kinase =~/0/){ #kinase gene
            if(exists $hash1{$symbol}){ #在network的gene
                print $O2 "$symbol\n";
                $hash2{$symbol}=1;
            }
            else{ #不在network的gene
                print $O1 "$symbol\n";
            }
        }
    }
}

close($O2);
my $line= readpipe("wc -l ../output/08_jintao_in_the_network.txt");
my @t= split/\s+/,$line;
my $true_line = $t[0];
my $the_number_of_jintao_in_network= $true_line -1;#减去header;
# print "$the_number_of_jintao_in_network\n";
# # print $true_line;

# # my @kinases;
while(<$I3>)
{
    chomp;
    my @f= split /\t/;
    unless(/^query/){
        my $query = $f[0];
        my $symbol = $f[1];
        if (exists $hash1{$symbol}){ #在network 的kinase
        print $O3 "$symbol\n";
        # print  "$symbol\n";
        my @neighbours=();
        my @tops;
        my $f4 ="../output/RWR/rwr_result_top/${symbol}.txt";
        open my $I4, '<', $f4 or die "$0 : failed to open input file '$f4' : $!\n";
            while(<$I4>)
            {
                chomp;
                my @f =split/\t/;
                my $neighbour = $f[0];
                push @tops,$neighbour;
                if (exists $hash2{$neighbour}){ #在top上的neighbour
                    push @neighbours,$neighbour;
                    # print "$symbol\t$neighbour\n";
                }  
            }
            # print "$symbol\t@neighbours\n";
            my $all_tops = join(";",@tops);
            my $all_events = join(";",@neighbours);
            print $O6 "$symbol\t$all_tops\n";
            print $O7 "$symbol\t$all_events\n";
            my $kinase_in_the_top = @neighbours;
            my $gene_not_kinase_in_the_top = 139 - $kinase_in_the_top; #TOP数目一共为139个。
            my $kinase_out_of_top = $the_number_of_jintao_in_network - $kinase_in_the_top;
            my $gene_not_kinase_out_of_the_top = 13852 -$gene_not_kinase_in_the_top - $kinase_out_of_top-$kinase_in_the_top;
            # my $kinase_in_the_top_final = $kinase_in_the_top -1; #因为top的第一个是自己，所以再减去一个
            my $output = "$symbol\t$kinase_in_the_top\t$gene_not_kinase_in_the_top\t$kinase_out_of_top\t$gene_not_kinase_out_of_the_top";
            print $O5 "$output\n";
        }
        else{ #不在network 的kinase
            print $O4 "$symbol\n";
        }
    }
}
