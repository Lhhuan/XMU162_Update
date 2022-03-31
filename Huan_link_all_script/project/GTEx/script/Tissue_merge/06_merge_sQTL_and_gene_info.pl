#为"/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz" 找到gene TSS并转bed得../output/Tissue_merge/Cis_eQTL/06_sig_eQTL_gene_TSS.bed.gz,排序得"../../output/Tissue_merge/Cis_eQTL/06_sQTL_gene_TSS.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;




my $fo1 = "../../output/Tissue_merge/Cis_eQTL/06_sQTL_gene_TSS.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Chr\tStart\tEnd\tegene\thotspot\ttissue\n";

my $f1 = "mart_export.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../../output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_sQTL_gene.txt.gz";
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    unless(/^Gene/){
        my $ENSG = $f[0];
        my $TSS =$f[1];
        my $start=$f[2];
        my $end=$f[3];
        my $chr =$f[4];
        my $TSS_S = $TSS -1; #1-based
        my $v = "chr${chr}\t$TSS_S\t$TSS";
        push @{$hash1{$ENSG}},$v;
    }
}

while(<$I2>)
{
    chomp;
    my @f =split/\t/;
    unless(/^SNP_chr/){
        my $SNP_chr = $f[0];
        my $SNP_pos =$f[1]; #0-based
        my $Pvalue=$f[2];
        my $gene_id=$f[3];
        $gene_id =~ s/\..*+//g;
        my $start = $SNP_pos;
        my $end =$SNP_pos +1;
        my $chr="chr${SNP_chr}";
        if($Pvalue <5E-8){
            # print "$f[3]\t$gene_id\n";
            if(exists $hash1{$gene_id}){
                my @vs=@{$hash1{$gene_id}};
                foreach my $v(@vs){
                    my $output= "$chr\t$start\t$end\t$gene_id\t$v";
                    unless(exists $hash2{$output}){
                        $hash2{$output}=1;
                        print $O1 "$output\n";
                    }
                }
            }
        }

    }
}

close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n| gzip >../../output/Tissue_merge/Cis_eQTL/06_sQTL_gene_TSS_sorted.bed.gz"