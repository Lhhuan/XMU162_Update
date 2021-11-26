#将../data/FIsInGene_122718_with_annotations.txt中的start和end用symbol表示（map到../output/03_final_all_network_gene_symbol.txt）
#得../output/04_network_symbol.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 ="../data/FIsInGene_122718_with_annotations.txt";
my $f2 ="../output/03_final_all_network_gene_symbol.txt";
my $fo1 ="../output/04_network_symbol.txt"; 
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

print $O1 "Symbol1\tSymbol2\tAnnotation\tDirection\tScore\n";
my(%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7);

while(<$I1>)
{
    chomp;
    my @f= split /\t/;
    unless(/^Gene1/){
        my $gene1 = $f[0];
        my $gene2 = $f[1];
        my $Annotation = $f[2];
        my $Direction = $f[3];
        my $Score = $f[4];
        my $v = "$gene2\t$Annotation\t$Direction\t$Score";
        push @{$hash1{$gene1}},$v;
    }
}

while(<$I2>)
{
    chomp;
    my @f= split /\t/;
     unless(/^Gene/){
         my $gene = $f[0];
         my $symbol = $f[1];
         $hash2{$gene} = $symbol;
     }
}

#对第一列基因进行id的转换
foreach my $ID (sort keys %hash1){
    if (exists $hash2{$ID}){
        my $s = $hash2{$ID};
        my @v =@{$hash1{$ID}};
        foreach my $v(@v){
            my $group1 = "$s\t$v";
            my @f= split/\t/,$group1;
            my $gene2 = $f[1];
            my $vg = join ("\t", @f[0,2,3,4]); #$f[0]是原来表格中的gene1转换后的id，也就是$s,
            push @{$hash3{$gene2}},$vg;
            #print STDERR "$group1";
        }
    }
}
# #对第二列基因进行转换
foreach my $ID (sort keys %hash3){
    if (exists $hash2{$ID}){
        my $s = $hash2{$ID};
        my @v =@{$hash3{$ID}};
        foreach my $v(@v){
            my $group2 = "$s\t$v";
            my @f=split/\t/,$group2;
            my $gene2 = $f[0];
            my $gene1 = $f[1];
            my $Annotation = $f[2];
            my $Direction = $f[3];
            my $Score = $f[4];
            print $O1 "$gene1\t$gene2\t$Annotation\t$Direction\t$Score\n";
        }
    }
}

            

