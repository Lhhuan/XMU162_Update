
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/eQTL_Catalogue/script/Knowledge_Graph/output/09_hotspot_5e_8egene_pos.bed.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/09_egene_pos_5e_8.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;
while(<$I1>)
{
    chomp;
    my @f=split/\t/;
    my $ensg =$f[3];
    my $chr=$f[4];
    my $start =$f[5];
    my $end =$f[6];
    my $output="$chr\t$start\t$end\t$ensg";
    unless(exists $hash1{$output}){
        $hash1{$output}=1;
        print $O1 "$output\n";
    }
    # if($chr=~/chr5/ & $start =~/181195/){
    #     print "$_\n";
    #     print "$output\n";
    # }
}

close($O1);
close($I1);
my $sorted_input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/09_egene_pos_5e_8_sorted.bed.gz"; 
system "zless $fo1 |sort -k1,1 -k2,2n |gzip>$sorted_input_file";
my $anno_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/";
my $out_dir = "/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/egene";


my $input_file_base_name = basename($sorted_input_file);

$ENV{'sorted_input_file'}  = $sorted_input_file; #设置环境变量
$ENV{'input_file_base_name'} = $input_file_base_name ;
$ENV{'output_dir'} = $out_dir ;
$ENV{'anno_dir'} = $anno_dir ;

system "bash annotation_marker_interval18_signalValue.sh";