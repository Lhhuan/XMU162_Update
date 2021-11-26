#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $command = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/non_factor/hotspot/0.01/promoter_whole_blood_segment_hotspot_cutoff_0.01.bed.gz | wc -l" ;


my $line_count = wc($command);
print "$line_count\n";

sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}




# my $median_rwr_normal_P_value = mid( @normal_score_P);

# sub mid{
#     my @list = sort{$a<=>$b} @_;
#     my $count = @list;
#     if( $count == 0 )
#     {
#         return undef;
#     }   
#     if(($count%2)==1){
#         return $list[int(($count-1)/2)];
#     }
#     elsif(($count%2)==0){
#         return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
#     }
# }  