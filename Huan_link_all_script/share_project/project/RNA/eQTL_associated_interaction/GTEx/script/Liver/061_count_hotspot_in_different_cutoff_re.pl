#统计不同cutoff下hotspot的数目/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18_re.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my @intervals = (18,15,12,9,8,7,6);
my @cutoffs = ();

for (my $i=0.01;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}


my $tissue = "Liver";
my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18_re.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Cutoff\tNumber_of_hotspot\n";

foreach my $cutoff(@cutoffs){
    my $hotspot_file = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my @arg1s = stat ($hotspot_file);
    my $hotspot_size = $arg1s[7]; #取$factor_hotspot_file size
    # print "$hotspot_size\n";
    # ------------------------------ factor_hotspot_size
    if ($hotspot_size >20){ #空compressed file is 20
        my $command_hotspot = "zless $hotspot_file | wc -l" ;
        my $hotspot_line_count = wc($command_hotspot);
        print $O1 "$cutoff\t$hotspot_line_count\n";
    }
    else{ #空 compressed file
        my $hotspot_line_count = 0;
        print $O1 "$cutoff\t$hotspot_line_count\n";
    }
}

sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}

