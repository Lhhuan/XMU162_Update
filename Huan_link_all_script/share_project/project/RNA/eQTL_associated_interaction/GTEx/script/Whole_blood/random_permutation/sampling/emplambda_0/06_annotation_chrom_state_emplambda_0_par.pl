#对10000 background_original_random/0进行chromatin_states的注释，得"/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2"/${i}_15_state_resemble_${input_file_base_name}
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
my $cutoff2 = "0";
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


# my $mark25 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/25_state_12_mark/${roadmap_biospecimen}_25_imputed12marks_mnemonics.bed.gz";
my $mark15 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/15_state_5_mark/${roadmap_biospecimen}_15_coreMarks_mnemonics.bed.gz";

my $anno_dir = "/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/background_original_random/$cutoff2";
#-----------------------------------------------------------------------------
my $input_dir=  "../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/background_original_random/${cutoff}/$cutoff2";
if(-e $anno_dir){
    print "${anno_dir}\texist\n";
}
else{
    system "mkdir -p $anno_dir";
}
my $pm = Parallel::ForkManager->new(40); ## 设置最大的线程数目
for (my $i=1;$i<10001;$i++){
# for (my $i=1;$i<101;$i++){
# for (my $i=500;$i<1001;$i++){
    my $pid = $pm->start and next; #开始多线程
    print "$i\n";
    my $input_file_base_name = "${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    my $sorted_input_file = "$input_dir/sorted_$input_file_base_name";

    my $command3 = "bedtools intersect   -a $sorted_input_file -b $mark15 -wo | gzip > $anno_dir/15_state_resemble_${input_file_base_name}";
    system $command3;

    $pm->finish;  #多线程结束
}
print "$type\tend\n";


