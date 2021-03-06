#对10000个../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2/*进行histone marker的注释，得/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2/*
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
my $type = "background_original_random";
my $cutoff2 = "0_0.176";
# my @types = ("random_permutation");
#----------------------- get roadmap_biospecimen
my $f1 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/tissue_id.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/GTEx_tissue/){
        my @f = split/\t/;
        my $GTEx_tissue = $f[0];
        my $Epigenomics_roadmap_biospecimen = $f[1];
        $hash1{$GTEx_tissue}=$Epigenomics_roadmap_biospecimen;
    }
}

my $roadmap_biospecimen = $hash1{$tissue};
print "$roadmap_biospecimen\n";

my $pm = Parallel::ForkManager->new(20); ## 设置最大的线程数目
# foreach my $type(@types){
my $out_dir = "/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2";
#-----------------------------------------------------------------------------
my $input_dir=  "../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2";

if(-e $out_dir){
    print "${out_dir}\texist\n";
}
else{
    system "mkdir -p $out_dir";
}
my $anno_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_biospecimen}";
for (my $i=1;$i<1001;$i++){
# for (my $i=1;$i<3;$i++){
    my $pid = $pm->start and next; #开始多线程
    my $input_file_base_name = "${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"; 
    my $input_file = "$input_dir/$input_file_base_name";
    my $sorted_input_file = "$input_dir/sorted_$input_file_base_name";
    my $command1 = "zless $input_file | sort -k1,1 -k2,2n | gzip > $sorted_input_file";
    $ENV{'sorted_input_file'}  = $sorted_input_file; #设置环境变量
    $ENV{'input_file_base_name'} = $input_file_base_name ;
    $ENV{'output_dir'} = $out_dir ;
    $ENV{'anno_dir'} = $anno_dir ;
    $ENV{'roadmap_biospecimen'} = $roadmap_biospecimen ;
    # #-----------annotation
    system $command1;
    my $command2 = "bash annotation_marker_interval18.sh";
    # my $command4 = "rm $input_file ";

    system $command2;
    # system $command4;

    print "$i\n";
    # print "$"
    $pm->finish;  #多线程结束
}
print "$type\tend\n";

# }

