#"../../output/${cancer}/Cis_eQTL/anno/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index.txt.gz";
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
# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");
# my @cancers =  ("KIRP"); 
my $f4 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/TCGA_mark.txt";
open my $I4, '<', $f4 or die "$0 : failed to open input file '$f4' : $!\n";  
my %hash1;
while(<$I4>)
{
    chomp;
    my @f = split/\t/;
    my $cancer = $f[0];
    my $marker = $f[1];
    push @{$hash1{$cancer}},$marker;
    push @{$hash1{$cancer}},"CHROMATIN_Accessibility";
    push @{$hash1{$cancer}},"TFBS";
}



my $anno_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell";
foreach my $cancer(@cancers){
    # my $f1 = "../../output/${cancer}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${cancer}.txt.gz";
    my $f1 = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $out_dir = "../../output/${cancer}/Cis_eQTL/anno";
    my $input_file_base_name = basename($f1);

    if(exists $hash1{$cancer}){
        my $fo1 = "$out_dir/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index.txt.gz";
        # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        open my $O1, "| gzip >$fo1" or die $!;

        print $O1 "hotspot_chr\thotspot_start\thotspot_end\tMarker\tjaacard_index\n";
        my @markers = @{$hash1{$cancer}};
        @markers=uniq(@markers); #数组内去重
        foreach my $marker(@markers){
            my $f2 = "${out_dir}/${marker}_${input_file_base_name}";
            if(-e $f2){ #判断 annotation文件存在
                my %hash2;
                # print "$f2\n";
                # print "$cancer\t$marker\n";
                open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                while(<$I2>)
                {
                    chomp;
                    my @f = split/\t/;
                    my $hotspot_chr = $f[0];
                    my $hotspot_start = $f[1];
                    my $hotspot_end = $f[2];
                    my $factor_chr = $f[3];
                    my $factor_start = $f[4];
                    my $factor_end = $f[5];
                    my $overlap_bp = $f[-1];
                    my $factor_length = $factor_end - $factor_start;
                    my $k = "$hotspot_chr\t$hotspot_start\t$hotspot_end";
                    my $v = "$factor_chr\t$factor_start\t$factor_end\t$overlap_bp";
                    push @{$hash2{$k}},$v;
                }
                
                my $hash_length = keys %hash2;
                # print "$hash_length\n";
                if($hash_length >0){ #annotation file is not null
                # print "$cancer\t$marker\n";
                    foreach my $k(sort keys %hash2){
                        my @vs = @{$hash2{$k}};
                        my %hash;
                        @vs = grep { ++$hash{$_} < 2 } @vs;
                        my @overlaps =();
                        my @factor_lengths =();
                        my @ss = split/\t/,$k;
                        my $hotspot_chr = $ss[0];
                        my $hotspot_start = $ss[1];
                        my $hotspot_end = $ss[2];
                        my $hotspot_length = $hotspot_end - $hotspot_start;
                        foreach my $v(@vs){
                            my @t=split/\t/,$v;
                            my $factor_chr = $t[0];
                            my $factor_start = $t[1];
                            my $factor_end =$t[2];
                            my $overlap_bp =$t[3];
                            my $factor_length = $factor_end - $factor_start;
                            push @overlaps, $overlap_bp;
                            push @factor_lengths, $factor_length;
                        }
                        my $all_overlap = sum @overlaps;
                        my $all_factor_length =sum @factor_lengths;
                        my $denominator = $all_factor_length+$hotspot_length-$all_overlap;
                        my $jaccard_index = $all_overlap/$denominator;
                        print $O1 "$k\t$marker\t$jaccard_index\n";
                    }
                    my $f3= $f1;
                    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

                    while(<$I3>)
                    {
                        chomp;
                        my @f = split/\t/;
                        my $chr = $f[0];
                        my $start = $f[1];
                        my $end =$f[2];
                        my $k = "$chr\t$start\t$end";
                        unless (exists $hash2{$k}){
                            print $O1 "$k\t$marker\t0\n";
                        }
                    }
                }
            }
        }
    }
}
