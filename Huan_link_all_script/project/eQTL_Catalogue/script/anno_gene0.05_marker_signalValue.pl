
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $sorted_input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene0.05_pos_sorted.bed.gz"; 
my $anno_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/";
my $out_dir = "/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/egene";


my $input_file_base_name = basename($sorted_input_file);

$ENV{'sorted_input_file'}  = $sorted_input_file; #设置环境变量
$ENV{'input_file_base_name'} = $input_file_base_name ;
$ENV{'output_dir'} = $out_dir ;
$ENV{'anno_dir'} = $anno_dir ;

system "bash annotation_marker_interval18_signalValue.sh";