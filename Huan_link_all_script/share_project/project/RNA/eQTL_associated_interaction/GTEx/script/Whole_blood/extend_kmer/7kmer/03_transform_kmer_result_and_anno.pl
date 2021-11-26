
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/7/kmer/extend/communities.csv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/7/kmer/extend/communities.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my %hash1;

while(<$I1>)
{
    chomp;
    if(/^>/){
        my @f = split/,/;
        my $p = $f[0];
        my $class =$f[1];
        $p=~s/>//g;
        # print "$p\n";
        my @ps= split/:/,$p;
        my $chr = $ps[0];
        my @ses = split/-/,$ps[1];
        my $start = $ses[0];
        my $end = $ses[1];
        print $O1 "$chr\t$start\t$end\t$class\n";

    }
}

close($I1);
close($O1);


$ENV{'input_file'}  = $fo1;
$ENV{'output_dir'}  = "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/7/kmer/extend/anno";

my $command = "bash annotation_factor_bedtools_intersect_interval18.sh";
system $command;
