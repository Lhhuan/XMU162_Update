#对10000个original_random和random_permutation进行chromatin_states的注释，得../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/${type}/${i}_25_state_resemble_${input_file_base_name}
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
my @types = ("random_permutation","original_random");
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
my $mark25 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/25_state_12_mark/${roadmap_biospecimen}_25_imputed12marks_mnemonics.bed.gz";
my $mark15 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/15_state_5_mark/${roadmap_biospecimen}_15_coreMarks_mnemonics.bed.gz";
foreach my $type(@types){
    my $anno_dir = "../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/${type}";
    #-----------------------------------------------------------------------------
    # my $input_file = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz"; #hotspot
    # my $input_file_base_name = basename($input_file);
    # $input_file_base_name =~ s/whole_blood/Whole_Blood/g;
    my $input_dir=  "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/${type}/${cutoff}";

    if(-e $anno_dir){
        print "${anno_dir}\texist\n";
    }
    else{
        system "mkdir -p $anno_dir";
    }

    my $genome="/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
    for (my $i=1;$i<10001;$i++){
    # for (my $i=500;$i<1001;$i++){
        my $pid = $pm->start and next; #开始多线程
        my $input_file_base_name = "${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
        my $input_file = "$input_dir/$input_file_base_name";
        my $sorted_input_file = "$input_dir/sorted_$input_file_base_name";
        my $command1 = "zless $input_file | sort -k1,1 -k2,2n | gzip > $sorted_input_file";
        # my $out_basename = 
        # print "$i\n";
        # #-----------annotation
        my $command2 = "bedtools intersect   -a $sorted_input_file -b $mark25 -wo | gzip > $anno_dir/25_state_resemble_${input_file_base_name}";
        my $command3 = "bedtools intersect   -a $sorted_input_file -b $mark15 -wo | gzip > $anno_dir/15_state_resemble_${input_file_base_name}";
        my $command4 = "rm $input_file ";
        system $command1;
        system $command2;
        system $command3;
        # system $command4;

        # print "$command1\n";
        # print "$command2\n";
        # print "$command3\n";
        # print "$command4\n";
        print "$i\n";
        $pm->finish;  #多线程结束
    }
    print "$type\tend\n";

}

