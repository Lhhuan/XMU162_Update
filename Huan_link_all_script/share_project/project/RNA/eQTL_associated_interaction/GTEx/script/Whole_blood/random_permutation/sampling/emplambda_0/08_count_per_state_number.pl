#对"$dir/filter_$name"中的chromatin state 进行统计，得"../../../output/${tissue}/Cis_eQTL/enrichment/interval_18/ALL/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Whole_Blood";


my $state = 15;
my $cutoff2 = "0";

my $output_dir = "../../../../../output/${tissue}/Cis_eQTL/enrichment/interval_18/ALL/" ;
my $fo1 = "$output_dir/${cutoff2}_${state}_state_count_${tissue}_cutoff_${cutoff}.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
print $O1 "chromatin_state\tnumber\trandom_number\n";
my $input_dir=  "../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2";
my $dir = "/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2";
for (my $i=1;$i<10001;$i++){  
    my $name = "${state}_state_resemble_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    my $f1 = "$dir/filter_$name";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件      
    my %hash1;

    while(<$I1>)
    {
        chomp;
        unless(/^hotspot_chr/){
        my @f =split/\t/;
        my $hotspot_chr=$f[0];
        my $hotspot_start = $f[1];
        my $hotspot_end = $f[2];
        my $chromatin_state = $f[6];
        my $v = join("\t",@f[0..2]);
        push @{$hash1{$chromatin_state}},$v;
        }
    }

    foreach my $k1 (sort keys %hash1){
    my @vs = @{$hash1{$k1}};
    my %hash2;
    my @uniq_vs = grep { ++$hash2{ $_ } < 2; } @vs;
    my $count = @uniq_vs;
    print $O1 "$k1\t$count\t$i\n";
    }
    print "$i\n";
}



