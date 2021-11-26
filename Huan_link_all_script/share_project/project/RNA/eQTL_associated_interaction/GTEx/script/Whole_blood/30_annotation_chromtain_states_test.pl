#interval_18 时，对../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz用annotation_chromatin_states_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
#extend 100bp
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my @intervals = (18,15,12,9,8,7,6);
# my @cutoffs = ();
# for (my $i=0.05;$i<0.96;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
#     push @cutoffs,$i;
#     # print "$i\n";
# }
# push @cutoffs,0.01;
# push @cutoffs,0.99;
# my @groups = ("hotspot","non_hotspot");
my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Whole_Blood";


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


my $input_file = "../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz";
# my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
my $input_file_base_name = basename($input_file);
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
# my $dir = dirname($script);
# print "$input_file_base_name\n";
my $output_dir = "../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}"; 
# mkdir $PMID;
# #------------
if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}
#------------
$ENV{'input_file'}  = $input_file; #设置环境变量
$ENV{'input_file_base_name'} = $input_file_base_name ;
$ENV{'output_dir'} = $output_dir ;
$ENV{'roadmap_biospecimen'} = $roadmap_biospecimen ;
# $ENV{'fraction'} = $fraction ;
my $command = "bash annotation_chromatin_states_interval18.sh";
# print "$command\n";
system $command;
print "$group\t$cutoff\n";




# my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
# foreach my $group(@groups){
#     foreach my $cutoff(@cutoffs){

#         my $pid = $pm->start and next; #开始多线程
        
#         $pm->finish;  #多线程结束
#     }
# } 


