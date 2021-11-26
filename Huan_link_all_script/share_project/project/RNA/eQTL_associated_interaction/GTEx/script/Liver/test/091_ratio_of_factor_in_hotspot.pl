##看各种factor在"~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}.txt.gz"中hotspot 在factor中的比例得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/091_factor_ratio_of_hotspot.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @cutoffs = ();
for (my $i=0.05;$i<0.96;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
push @cutoffs,0.01;
# push @cutoffs,0.99;

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/091_factor_ratio_of_hotspot.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Cutoff\tFactor\tNumber_of_hotspot_in_factor\tNumber_of_hotspot\tFactor_ratio\n";  
foreach my $cutoff(@cutoffs){
    my $f1 = "~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    #--------------------
    my (%hash1,%hash2,%hash3,%hash4);
    #-------hotspot
    # print "$f1\n";
    my @factors=();
    my @hotspots= ();
    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $chr = $f[0];
        my $start = $f[1];
        my $end = $f[2];
        # my @factors = @f[3..12];
        my $k = join("\t",@f[0..2]);
        push @hotspots,$k;
        if(/\bchr\b/){
        #    print $O1 "chr\tstart\tend\tscore\n";
            @factors = @f;
        }
        else{
            for (my $i=3;$i<5;$i++){
                if (abs($f[$i] -0) < 0.000000001){ # hotspot:0
                    push @{$hash1{$factors[$i]}},$f[$i];
                }
                else{
                    push @{$hash2{$factors[$i]}},$f[$i];# hotspot:1
                }
            }
        }
    }

    my $number_of_hotspot = @hotspots;
    foreach my $factor(sort keys %hash1){
        my @v1s = @{$hash1{$factor}};
        my @n2s=();
        if (exists $hash2{$factor}){ # hotspot:1 TP
            my @vs2 = @{$hash2{$factor}};
            my $number_of_hotspot_in_factor =@vs2;
            my $factor_ratio= $number_of_hotspot_in_factor/$number_of_hotspot;
            $factor_ratio=sprintf "%.4f",$factor_ratio;
            print $O1 "$cutoff\t$factor\t$number_of_hotspot_in_factor\t$number_of_hotspot\t$factor_ratio\n";
        }
        else{
            my $number_of_hotspot_in_factor =0;
            my $factor_ratio= $number_of_hotspot_in_factor/$number_of_hotspot;
            $factor_ratio=sprintf "%.4f",$factor_ratio;
            print $O1 "$cutoff\t$factor\t$number_of_hotspot_in_factor\t$number_of_hotspot\t$factor_ratio\n";
        }
    }


}


