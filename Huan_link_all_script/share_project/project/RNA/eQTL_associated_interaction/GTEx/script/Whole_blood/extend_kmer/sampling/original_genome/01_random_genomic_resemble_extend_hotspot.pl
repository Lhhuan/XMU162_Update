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
my $group ="hotspot";
my $tissue = "Whole_Blood";
my $cutoff2 = "original_genome";

my $input_file = "../../../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6/extend_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
my $input_file_base_name = basename($input_file);
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
my $output_dir=  "/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/$cutoff2/bed";
my $fa_out_dir ="/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/$cutoff2/fa";

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
    my $out_file_base_name = "${i}_resemble_${input_file_base_name}";
    $out_file_base_name =~ s/bed\.gz/fa/;
    
    my $command2 = " bedtools getfasta -fi ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa -bed  $out_file -fo $fa_out_dir/$out_file_base_name";
    system $command2;
    # print "$i\n";
    print "$i\n";
}

