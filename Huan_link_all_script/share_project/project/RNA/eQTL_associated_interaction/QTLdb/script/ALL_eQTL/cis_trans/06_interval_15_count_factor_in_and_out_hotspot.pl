#看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_${group}_${type}.txt.gz"中hotspot 和 non-hotspotfactor的个数，得
#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/factors_score/06_annotation_merge_factors_score_${group}_${type}.txt.gz"
#0-based
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my @types= ("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
my @groups = ("hotspot","non_hotspot");


foreach my $type(@types){
    foreach my $group(@groups){
        my $f1 = "~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_${group}_${type}.txt.gz";
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        #--------------------
        my $fo1 = "~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/factors_score/06_annotation_merge_factors_score_${group}_${type}.txt.gz";
        # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        open my $O1, "| gzip >$fo1" or die $!;
        while(<$I1>)
        {
            chomp;
            my @f = split/\t/;
            my $chr = $f[0];
            my $start = $f[1];
            my $end = $f[2];
            my @factors = @f[3..12];
            
            if(/\bchr\b/){
               print $O1 "chr\tstart\tend\tscore\n";
            }
            else{
                my $score = sum @factors;
                # print "@factors\n";
                print $O1 "$chr\t$start\t$end\t$score\n";
            }
        }
    }

}
