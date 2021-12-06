
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "./output/03_qtl_p_emplamda.bed.gz"; 
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;

my %hash1;


while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $k=join("\t",@f[0,1]);#pos,chr
    my $p = $f[2];
    push @{$hash1{$k}},$p;
}
print "read2\n";

while(<$I2>)
{
    chomp;
    unless(/^emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $chr =$f[2];
        my $pos= $f[1];
        my $end= $pos+1;
        my $start = $end-1;
        my $k="$chr\t$start";
        if(exists $hash1{$k}){
            $chr="chr${chr}";
            my @vs= $hash1{$k};
            my @sorted_vs = sort {$a <=> $b} @vs;
            my $v= $sorted_vs[0]; #min p value
            if($v< 5E-8){
                print $O1 "$chr\t$start\t$end\t$v\tT\t$emplambda\n";
            }
            else{
                print $O1 "$chr\t$start\t$end\t$v\tF\t$emplambda\n";
            }
            

        }
        else{
            print "$k\n";
        }

    }
}

close($O1);
# system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ./output/03_hotspot_without_egene_sorted.bed.gz";




