
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;


my $anno_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/";
my $out_dir = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/annotation/extend";

my $sorted_input_file = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz";
my $input_file_base_name = "Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz";

$ENV{'sorted_input_file'}  = $sorted_input_file; #设置环境变量
$ENV{'input_file_base_name'} = $input_file_base_name ;
$ENV{'output_dir'} = $out_dir ;
$ENV{'anno_dir'} = $anno_dir ;

system "bash annotation_marker_interval18.sh";