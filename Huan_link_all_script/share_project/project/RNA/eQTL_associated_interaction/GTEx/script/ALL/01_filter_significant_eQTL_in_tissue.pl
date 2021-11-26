#过滤"../../data/GTEx_Analysis_v8_eQTL_hg19/${tissue}${suffix}"中显著的eQTL，得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/significant_eQTL.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager;

# my @tissues = ("Whole_Blood","Lung");
my @tissues = ("Whole_Blood","Lung","Ovary","Liver","Pancreas");

foreach my $tissue(@tissues){
# my $tissue = "Whole_Blood"; 
    my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";

    my $f1 = "../../data/GTEx_Analysis_v8_eQTL_hg19/${tissue}${suffix}";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my(%hash1,%hash2,%hash3,%hash4);

    my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/significant_eQTL.txt.gz";
    open my $O1, "| gzip >$fo1" or die $!;

    while(<$I1>)
    {
        chomp;
        unless(/^variant_id/){
            my @f = split/\t/;
            my $SNP_chr =$f[1];
            my $SNP_pos =$f[2];
            my $Pvalue =$f[-6];
            my $start = $SNP_pos;
            my $end = $SNP_pos+1;
            if ($Pvalue <5e-8){
                my $output = "chr$SNP_chr\t$SNP_pos\t$end";
                # print "$output\n";
                unless (exists $hash1{$output}){
                    $hash1{$output} =1;
                    print $O1 "$output\n";
                }
            }
        }
    }
}

