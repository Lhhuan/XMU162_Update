#产生1000个与../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}",进行chromatin  annotation 得$anno_dir/${i}_25_state_resemble_${input_file_base_name}
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

my $anno_dir = "../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/original_random";
my $mark25 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/25_state_12_mark/${roadmap_biospecimen}_25_imputed12marks_mnemonics.bed.gz";
my $mark15 = "/share/data0/GTEx/annotation/ROADMAP/chromatin_states/15_state_5_mark/${roadmap_biospecimen}_15_coreMarks_mnemonics.bed.gz";
#-----------------------------------------------------------------------------
my $input_file = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz"; #hotspot
my $input_file_base_name = basename($input_file);
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
my $output_dir=  "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/original_random/${cutoff}";

if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}

my $genome="/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
for (my $i=1;$i<1001;$i++){
    my $out_file = "$output_dir/${i}_resemble_${input_file_base_name}";
    #generate random file 
    my $command1 = "bedtools shuffle -i $input_file -g $genome -excl $input_file -chrom | gzip >$out_file"; 
    # print "$command\n";
    system $command1;
    print "$i\n";
    # #-----------annotation
    my $command2 = "bedtools intersect   -a $out_file -b $mark25 -wo | gzip > $anno_dir/${i}_25_state_resemble_${input_file_base_name}";
    my $command3 = "bedtools intersect   -a $out_file -b $mark15 -wo | gzip > $anno_dir/${i}_15_state_resemble_${input_file_base_name}";
    system $command2;
    system $command3;
}



