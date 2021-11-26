#将../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 按照特定距离将emplambda取平均
#得../output/average_emplambda_xQTL.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;




my $f1 = "test.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "test.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

#---------------------------------
# my $f1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz";
# # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

# my $f2 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz";
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
#------------------------------------

my $fo1 = "../output/average_emplambda_xQTL.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "chr\tstart\tend\taverage_emplambda\txQTL\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $t =$f[1]; 
        my $pos = $t; #snp pos
        my $chr =$f[2]; #snp chr
        my $xQTL = $f[3];
        my $k1 = $chr;
        unless ($emplambda =~/NA/){ 
            push @{$hash1{$chr}},$pos;
            $hash2{$xQTL}=1;
        }
    }
}

my %hash3;
foreach my $k1(sort keys %hash1){
    my @pos = @{$hash1{$k1}};
    my $max_pos = max @pos;
    my $min_pos =min @pos;
    my $v = "$min_pos\t$max_pos";
    $hash3{$k1}=$v;
}

foreach my $chr(sort keys %hash3){
    my $v = $hash3{$chr};
    my @ts =split/\t/,$v;
    my $min_pos = $ts[0];
    my $max_pos = $ts[1];
    # print "$chr\n";

    foreach my $xQTL(sort keys %hash2 ){
        my $min_up = $max_pos-1000; #min 的上限
        for(my $i = $min_pos;$i<=$min_up; $i = $i+1000){
            my $min =$i;
            my $max = $i +1000;
            my @pos_interval =();
            print "$chr\n";
            # my $f2 = "test.txt";
            # open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
            while(<$I2>)
            {
                chomp;
                unless(/emplambda/){
                    my @f = split/\t/;
                    my $emplambda =$f[0];
                    my $t =$f[1]; 
                    my $pos = $t; #snp pos
                    my $chr_f =$f[2]; #snp chr
                    # print "$chr_f\n";
                    my $xQTL_f = $f[3];
                    # print "$chr\n";
                    if($chr_f =~ /\b$chr\b/){
                        # print "$chr_f\n";
                    # if($chr_f == $chr ){
                        # print "$chr_f\n";
    #                     # if ($xQTL_f =~/$xQTL/ ){
    #                     #     print "$chr\t$chr_f\t$xQTL\t$xQTL_f\n";
    #                     #     if ($pos >= $min && $pos < $max){
    #                     #         push @pos_interval, $emplambda;
    #                     #         # print "$emplambda\n";
    #                     #     }
    #                     # }
                    }
                
                }
            }
            # my $length = @pos_interval;
            # print "$length\t@pos_interval\n";
            # if ($length >0){
            #     my $all_emplambda = sum @pos_interval;
            #     my $average_emplambda = $all_emplambda/$length;
            #     print $01 "$chr\t$min\t$max\t$average_emplambda\t$xQTL\n";

            # }
            # else{
            #     my $average_emplambda =0;
            #     print $01 "$chr\t$min\t$max\t$average_emplambda\t$xQTL\n";
            # }
        }
        
    }
}
