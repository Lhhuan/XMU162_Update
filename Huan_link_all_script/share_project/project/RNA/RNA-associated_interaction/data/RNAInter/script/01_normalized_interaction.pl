#normalize ../output/${geneType}.txt, 得../normalized/01_${geneType}_normalized.txt;
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my %hash2;
my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";

my @genetypes = ("RNA-Protein","RNA-RNA","RNA-DNA");

foreach my $geneType(@genetypes){
    my $f2 = "../output/${geneType}.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
    print STDERR "$f2\n";
    my $fo1 = "../normalized/01_${geneType}_normalized.txt";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
    print $O1 "$header\n";
    while(<$I2>) 
    {
        chomp;
        unless(/^RNAInter/){
            my @f= split/\t/;
            my $RNAInter_ID = $f[0];
            my $Interactor1 = $f[1];
            my $ID1 =$f[2];
            my $Category1 =$f[3];
            my $Species1 =$f[4];
            my $Interactor2 =$f[5];
            my $ID2 = $f[6];
            my $Category2 =$f[7];
            my $Species2 =$f[8];
            my $Score =$f[9];
            if($Species1 =~/Homo sapiens/ && $Species2 =~/Homo sapiens/ ){
                my $chr1 = "NA";
                my $start1 = "NA";
                my $end1 = "NA";
                my $Strand1 = "NA";
                my $Type1 =$Category1;
                my $interaction_name1 = "$Interactor1";
                my $geneID1 = $ID1;
                my $geneName1 = "NA";
                my $chr2 = "NA";
                my $start2 ="NA";
                my $end2 ="NA";
                my $Strand2 = "NA";
                my $Type2 = $Category2;
                my $interaction_name2 =$Interactor2;
                my $GeneID2 = $ID2;
                my $geneName2 = "NA";
                my $Sub_interaction_type = "${Type1}-${Type2}";
                my $output = "$chr1\t$start1\t$end1\t$Strand1\t$Type1\t$interaction_name1\t$geneID1\t$geneName1\tNA\t$chr2\t$start2\t$end2\t$Strand2\t$Type2\t$interaction_name2\t$GeneID2\t$geneName2\tNA\t$Sub_interaction_type\t$geneType\tNA\tNA\tNA\tRNAInter";
                unless (exists $hash2{$output}){
                    $hash2{$output} =1;
                    print $O1 "$output\n";
                }
            }
        }
    }
}
