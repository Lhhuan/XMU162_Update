#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
# use List::Util qw/max min/;
# use Parallel::ForkManager;


my $f1 = "/home/qinwei/project/1st_DNAseq/vcf/annovar_anno/config";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "final_result.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;



my @vs = ();
while(<$I1>)
{
    chomp;
   my $f2_name= $_;
   push @vs,$f2_name;
}

my %hash1;

foreach my $f2_name(@vs){
   my $f2 = "/home/qinwei/project/1st_DNAseq/vcf/annovar_anno/${f2_name}";
#    print "$f2\n";
   open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        if(/^Chr/){
            my $output="$_\tClassification\tDP\tDP4\tGenotype\tAD\tVAF\tSample";
            unless(exists $hash1{$output}){
                # (exists $hash1{$k1})
                $hash1{$output} =1;
                print $O1 "$output\n";
            }
        }
        else{
            my $org_output = join("\t",@f[0..20]);
            my $data= $f[-1];
            my $describition =$f[-2];
            my $type_info =$f[-3];
            my @t =split/\:/,$data;
            my $Genotype =$t[0];
            my $AD = $t[2];
            # print "$AD\n";
            my @aa = split/\,/,$AD;
            my $ref= $aa[0];
            my $alt = $aa[1];
            my $VAF = $alt/($ref+$alt);
            my @dd = split/\;/,$type_info;
            my @DPs=();
            my @DP4s=();
            foreach my $d (@dd){
                if ($d =~/DP=/){
                    push @DPs,$d;
                }
                elsif($d =~/DP4=/){
                    push @DP4s,$d;
                }
            }
            my $DP = $DPs[0];
            my $DP4 = $DP4s[0];
            $DP =~ s/DP=//g;
            $DP4 =~ s/DP4=//g;
            my $sample =$f2_name;
            $sample =~ s/\.hg19_multianno.txt//g;
            my $output1 = $org_output;
            my $output2 =  "$DP\t$DP4\t$Genotype\t$AD\t$VAF\t$sample";
            if ($type_info =~/INDEL/){
                print $O1 "$output1\tINDEL\t$output2\n";
            }
            else{
                print $O1 "$output1\tSNP\t$output2\n";
            }
        }
    }
}