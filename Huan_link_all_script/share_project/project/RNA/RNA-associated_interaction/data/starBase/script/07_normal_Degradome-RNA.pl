#normalize ../raw_data/ENCORI_hg19_degradome-seq_all_Degradome-${geneType}_Degradome-RNA.txt, 得../normalized/07_ENCORI_hg19_degradome-seq_all_normalized.txt,
#得normalized fail的文件../output/07_ENCORI_hg19_degradome-seq_all_normalized_fail.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../../miRBase/01_hsa_hg19.gff3";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../normalized/07_ENCORI_hg19_degradome-seq_all_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../output/07_ENCORI_hg19_degradome-seq_all_normalized_fail.txt"; #得没有匹配的mirna文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my %hash1;
my %hash2;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless (/^Chr/){
        my $chr = $f[0];
        my $pos1 = $f[3];
        # my $start = $pos1 -1; #gff文件为1 based, 之后的文件用0-based
        my $start = $pos1;
        my $pos2 = $f[4];
        my $end = $pos2;
        my $Strand = $f[6];
        my $Info = $f[8];
        my @ts= split/;/,$Info;
        my $miRNA_ID= $ts[0];
        $miRNA_ID =~ s/ID=//g;
        foreach my $t(@ts){
            if ($t =~ /^Name/){
                $t =~ s/Name=//g;
                # print STDERR "$t\n";
                my $k = "$t\t$miRNA_ID";
                my $v = "$chr\t$start\t$end\t$Strand\tmiRNA\t$t\tNA\tNA\tNA";
                push @{$hash1{$t}},$v;
            }
        }
    }
}

my @genetypes = ("mRNA","ncRNA");
foreach my $geneType(@genetypes){
    my $f2 = "../raw_data/ENCORI_hg19_degradome-seq_all_Degradome-${geneType}_Degradome-RNA.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
    while(<$I2>) #用end - start确定该文件为0-based
    {
        chomp;
        unless(/^#/){
            my @f= split/\t/;
            if (/^miRNAid/){
                print $O2 "$_\n";
            }
            else{
                my $miRNAid =$f[0];
                my $miRNAname = $f[1];
                my $k = "$miRNAname\t$miRNAid";
                my $ENSG_ID = "NA";
                my $gene_symbol = $f[2];
                my $gene_type= $f[3];
                my $chr = "NA";
                # my $narrowStart ="NA";
                # my $narrowEnd ="NA";
                my $broadStart ="NA";
                my $broadEnd = "NA";
                my $strand ="NA";
                my $pancancerNum ="NA";
                my $tub_interaction_type = "miRNA-${geneType}";
                my $interaction_name2 = $gene_symbol;
                my $tissue = "NA";
                if (exists $hash1{$k}){ #
                    my @vs =@{$hash1{$k}};
                    foreach my $v(@vs){
                        my $output = "$v\t$chr\t$broadStart\t$broadEnd\t$strand\t${geneType}\t$interaction_name2\t$ENSG_ID\t$gene_symbol\tNA\t$tub_interaction_type\tRNA-RNA\tHg19\t$tissue\t$pancancerNum\tstarBase";
                        unless(exists $hash2{$output}){
                            $hash2{$output}=1;
                            print $O1 "$output\n";
                        }
                        
                    }
                }
                else{
                    my $v = "NA\tNA\tNA\tNA\tmiRNA\t$miRNAname\tNA\tNA\tNA";
                    my $output = "$v\t$chr\t$broadStart\t$broadEnd\t$strand\t${geneType}\t$interaction_name2\t$ENSG_ID\t$gene_symbol\tNA\t$tub_interaction_type\tRNA-RNA\tHg19\t$tissue\t$pancancerNum\tstarBase";
                    unless(exists $hash2{$output}){
                        $hash2{$output}=1;
                        print $O1 "$output\n";
                    }
                    print $O2 "$_\n";
                }

            }
        }
    }
}