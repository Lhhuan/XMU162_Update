#得"../../output/cancer_total/number_of_cancer_all_hotspot.txt.gz"，"../../output/cancer_total/number_of_cancer_specific_hotspot.txt.gz"
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


 
my $fo1 = "../../output/cancer_total/number_of_cancer_all_hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo3 = "../../output/cancer_total/number_of_cancer_specific_hotspot.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;

print $O3 "cancer\tNumber\n";
print $O1 "cancer\tNumber\n";
 
foreach my $cancer(@cancers){
    my $all_hotspot = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $cancer_specific_hotspot = "../../output/cancer_total/specific/pure/${cancer}/${cancer}_specific.bed.gz";

    my $command_all_hotspot = "zless $all_hotspot | wc -l" ;
    my $all_hotspot_line_count = wc($command_all_hotspot);
    print $O1 "$cancer\t$all_hotspot_line_count\n";

    my $command_cancer_specific_hotspot = "zless $cancer_specific_hotspot  | wc -l" ;
    my $cancer_specific_hotspot_line_count = wc($command_cancer_specific_hotspot);
    print $O3 "$cancer\t$cancer_specific_hotspot_line_count\n";
    print "$cancer\n";
}



sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
