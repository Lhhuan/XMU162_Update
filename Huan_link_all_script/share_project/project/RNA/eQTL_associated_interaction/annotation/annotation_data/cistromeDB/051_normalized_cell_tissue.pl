

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;
# use main::uniq;

my $f1 = "./normal_cell/05_normal_chrom_marker_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo4 = "./normal_cell/051_normal_chrom_marker_info_refine.txt";
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";

print $O4 "Marker\tDCid\tSpecies\tGSMID\tFactor\tCell_type\tTissue_type\tCell_line\tRefine_Cell_type\tRefine_tissue_type\n";

my $fo5 = "./normal_cell/051_normal_chrom_marker_info_refine_unique_cell.txt";
open my $O5, '>', $fo5 or die "$0 : failed to open output file '$fo5' : $!\n";
print $O5 "Cell_line\tRefine_Cell_type\tRefine_tissue_type\n";


my (%hash1,%hash2,%hash3,%hash4,%hash5);
while(<$I1>)
{
    chomp;
    unless(/^Marker/){
        my @f = split/\t/;
        my $Marker = $f[0];
        my $DCid =$f[1];
        my $Species =$f[2];
        my $GSMID =$f[3];
        my $Factor =$f[4];
        my $Cell_line =$f[5];
        my $Cell_type =$f[6];
        my $Tissue_type =$f[7];
        push @{$hash1{$Cell_line}},$Cell_type;
        push @{$hash2{$Cell_line}},$Tissue_type;
        my $v5 = "$Marker\t$DCid\t$Species\t$GSMID\t$Factor\t$Cell_type\t$Tissue_type";
        push @{$hash5{$Cell_line}},$v5;
    }
}

foreach my $k(sort keys %hash1){
    my @v1s = @{$hash1{$k}};
    my @v2s = @{$hash2{$k}};
    # my @uniq_v1 =uniq(@v1s);
    # my @uniq_v2 = uniq(@v2s);
    my (%ha1,%ha2);
    my @uniq_v1 = grep{++$ha1{$_}<2}@v1s;
    my @uniq_v2 = grep{++$ha2{$_}<2}@v2s;
    my $Cell_type = $uniq_v1[0];
    my $num2 = @uniq_v2;
    unless($k=~/\bV6.5\b/){ #	Mus musculus cell line
        if($k=~/\bAB32\b/){
            my $v3 = "$Cell_type\tBreast";
            $hash3{$k}=$v3;  
        }
        elsif($k=~/\bAG18359\b/){
            my $v3 = "$Cell_type\tBlood";
            $hash3{$k}=$v3;          
        }
        elsif($k=~/\bBG01V - ATCC SCRC-2002|KhES-1|RUES2|SA121|SEES-3|WA01|WIS2|hESC|VAL-3\b/){
            my $v3 = "$Cell_type\tEmbryo";
            $hash3{$k}=$v3;  
        }
        else{
            if($num2 <2){
                my $Tissue_type =$uniq_v2[0];
                if($Cell_type =~/iPSC/){
                    my $v3 = "$Cell_type\tiPSC";
                    $hash3{$k}=$v3;            
                }
                elsif($Cell_type =~/Lymphoblastoid/){
                    my $v3 = "$Cell_type\tBlood";
                    $hash3{$k}=$v3;              
                }
                else{
                    my $v3 = "$Cell_type\t$Tissue_type";
                    $hash3{$k}=$v3;
                }
            }
            else{
                foreach my $v2(@uniq_v2){
                    unless($v2 =~/None/){
                        my $Tissue_type =$uniq_v2[0];
                        my $v3 = "$Cell_type\t$Tissue_type";
                        $hash3{$k}=$v3; 
                    }
                }
            }
        }
    }
}


foreach my $k(sort keys %hash3){
    my $cell_tissue = $hash3{$k};
    print $O5 "$k\t$cell_tissue\n";
    my @infos = @{$hash5{$k}};
    foreach my $info(@infos){
        print $O4 "$info\t$k\t$cell_tissue\n";
    }
}