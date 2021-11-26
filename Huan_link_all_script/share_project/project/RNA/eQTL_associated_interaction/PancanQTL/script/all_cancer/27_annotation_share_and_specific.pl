#"../../output/${cancer}/Cis_eQTL/anno/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index_lable.txt.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;
my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 

my %hash1;
$hash1{"ACC"}="Adrenal_Gland";
$hash1{"BRCA"}="Breast_Mammary_Tissue";
#------------
$hash1{"KIRC"}="Kidney_Cortex";
$hash1{"KIRP"}="Kidney_Cortex";
$hash1{"LAML"}="Whole_Blood";
$hash1{"LIHC"}="Liver";
$hash1{"LUAD"}="Lung";
$hash1{"LUSC"}="Lung";
$hash1{"OV"}="Ovary";

#---------------------

$hash1{"PAAD"}="Pancreas";
$hash1{"PRAD"}="Prostate";
$hash1{"STAD"}="Stomach";
$hash1{"THCA"}="Thyroid";

foreach my $cancer(@cancers){
    my $out_dir = "../../output/${cancer}/Cis_eQTL/anno";
    my $f2 = "$out_dir/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index.txt.gz";
    if(-e $f2){
        if (exists $hash1{$cancer}){
            my %hash2;
            my $tissue = $hash1{$cancer};
            my $specific_file= "../../output/cancer_total/specific/pure/${cancer}/cancer_${cancer}_${tissue}_specific.bed.gz";
            my $f1 = $specific_file;
            open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
            my $fo1 = "$out_dir/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index_label.txt.gz";
            open my $O1, "| gzip >$fo1" or die $!;
            print $O1 "Chr\tstart\tend\tmarker\tjaccard_index\tClass\n";
            while(<$I1>)
            {
                chomp;
                my @f = split/\t/;
                my $chr =$f[0];
                my $start = $f[1];
                my $end =$f[2];
                my $k= join("\t",@f);
                $hash2{$k}=1;
            }
            #-------------------
            #--------------
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            while(<$I2>)
            {
                chomp;
                unless(/^hotspot_chr/){
                    my @f = split/\t/;
                    my $chr =$f[0];
                    my $start = $f[1];
                    my $end =$f[2];
                    my $k= join("\t",@f[0..2]);
                    my $marker = $f[3];
                    my $jaccard_index = $f[4];
                    if (exists $hash2{$k}){
                        print $O1 "$_\tSpecific\n";
                    }
                    else{
                        print $O1 "$_\tOthers\n";
                    }
                }
            }
        }
    }
}

