#得"../../output/cancer_total/08_number_of_cancer_all_hotspot.txt.gz"，"../../output/cancer_total/08_number_of_cancer_specific_hotspot.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $cancer = "Lung";
my $j = 18;

my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 


 
my $fo1 = "../../output/cancer_total/08_number_of_cancer_all_hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo3 = "../../output/cancer_total/08_number_of_cancer_specific_hotspot.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;

print $O1 "cancer\tNumber\n";
print $O3 "cancer\ttissue\tNumber\n";
 
foreach my $cancer(@cancers){
    my $all_hotspot = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $command_all_hotspot = "zless $all_hotspot | wc -l" ;
    my $all_hotspot_line_count = wc($command_all_hotspot);
    print $O1 "$cancer\t$all_hotspot_line_count\n";
}

my $f1 = "../../output/cancer_total/specific/pure/cancer_specific.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;

while(<$I1>)
{
    chomp;
    unless(/^h_chr/){
       my @f= split/\t/;
       my $v =join("\t",@f[0..2]);
       my $k = join("\t",@f[3,4]);
       push @{$hash1{$k}},$v;
    }
}    

foreach my $k(sort keys %hash1 ){
    my @vs = @{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs; 
    my $count =@vs;
    print $O3 "$k\t$count\n";
}







sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
