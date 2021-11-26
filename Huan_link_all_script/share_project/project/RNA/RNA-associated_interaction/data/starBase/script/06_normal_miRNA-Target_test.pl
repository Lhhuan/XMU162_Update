#normalize "../raw_data/ENCORI_hg19_CLIP-seq_all_miRNA-${geneType}_${program}_miRNA-Target.txt", 得../normalized/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized.txt,
#得normalized fail的文件../output/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized_fail.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../../miRBase/01_hsa_hg19.gff3";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../normalized/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../output/06_ENCORI_hg19_CLIP-seq_all_miRNA_target_normalized_fail.txt"; #得没有匹配的mirna文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";

my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";
my (%hash1,%hash2,%hash3);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless (/^Chr/){
        my $chr = $f[0];
        my $pos1 = $f[3];
        # my $start = $pos1 -1; #gff文件为1 based, 之后的文件用1-based
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
                my $k = "$t\t$miRNA_ID";
                my $v = "$chr\t$start\t$end\t$Strand\tmiRNA\t$t\tNA\tNA\tNA";
                # print STDERR "$t\n";
                push @{$hash1{$k}},$v;
            }
        }
    }
}

my @genetypes = ("mRNA", "lncRNA", "pseudogene", "circRNA", "sncRNA");
my @programs =("PITA,RNA22,miRmap","DIANA-microT,miRanda","PicTar,TargetScan");
foreach my $geneType(@genetypes){
    foreach my $program(@programs){
        my $f2 = "../raw_data/ENCORI_hg19_CLIP-seq_all_miRNA-${geneType}_${program}_miRNA-Target.txt";
        open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
        while(<$I2>) #用end - start确定该文件为1-based
        {
            chomp;
            unless(/^#/){
                my @f= split/\t/;
                my $number = @f;
                if ($number >1){
                    if (/^miRNAid/){
                        print $O2 "$_\tSub_interaction_type\n";
                    }
                    else{
                        my $miRNAid =$f[0];
                        my $miRNAname = $f[1];
                        my $k = "$miRNAname\t$miRNAid";
                        my $ENSG_ID = $f[2];
                        my $gene_name = $f[3];
                        # my $geneType =$f[4];
                        my $chr = $f[5];
                        my $narrowStart =$f[6];
                        my $narrowEnd =$f[7];
                        my $broadStart =$f[8];
                        my $broadEnd = $f[9];
                        my $strand =$f[10];
                        my $pancancerNum =$f[-1];
                        my $Sub_interaction_type = "miRNA-${geneType}";
                        my $interaction_name2 = $gene_name;
                        my $tissue = "NA";
                        if (exists $hash1{$k}){ #
                            my @vs =@{$hash1{$k}};
                            foreach my $v(@vs){
                            # if ($pancancerNum>1){ #pancancerNum 为minimum number of Cancer types 
                                
                                my $output = "$v\t$chr\t$broadStart\t$broadEnd\t$strand\t$geneType\t$interaction_name2\t$ENSG_ID\t$gene_name\tNA\t$Sub_interaction_type\tRNA-RNA\tHg19\t$tissue\t$pancancerNum\tstarBase";
                                unless(exists $hash2{$output}){
                                    $hash2{$output} =1;
                                    print $O1 "$output\n";
                                }
                            }
                        }
                        else{
                            my $v = "NA\tNA\tNA\tNA\tmiRNA\t$miRNAname\tNA\tNA\tNA";
                            my $output = "$v\t$chr\t$broadStart\t$broadEnd\t$strand\t$geneType\tinteraction_name2\t$ENSG_ID\t$gene_name\tNA\t$Sub_interaction_type\tRNA_RNA\tHg19\t$tissue\t$pancancerNum\tstarBase";
                            unless(exists $hash2{$output}){
                                $hash2{$output} =1;
                                print $O1 "$output\n";
                            }
                            my $output2 = "$_\t$Sub_interaction_type";
                            unless(exists $hash3{$output2}){
                                print $O2 "$output2\n";
                            }
                        }
                    }
                }
            }
        }
    }
}
