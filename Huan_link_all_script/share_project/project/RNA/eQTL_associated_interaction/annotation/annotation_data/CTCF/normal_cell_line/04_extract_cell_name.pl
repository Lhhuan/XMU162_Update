#提取../data/download.sh 中的cell line信息，得./unique_cell_line.txt
#利用"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_existing_ccle_cell_line_info.txt"， "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt";和"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt"提供cell line 的disease,得有disease ./04_cell_line_in_cistrome.txt,得没有disease文件 ./04_cell_line_without_cistrome.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_existing_ccle_cell_line_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $f3 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
my $f4 = "../data/download.sh";
open my $I4, '<', $f4 or die "$0 : failed to open input file '$f4' : $!\n"; 

# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#-------------------
my $fo1 = "./04_unique_cell_line.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $fo2 = "./04_cell_line_in_cistrome.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

my $fo3 = "./04_cell_line_without_cistrome.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";



print $O2 "Cell_line\tdisease\n";

my %hash1;
my %hash2;

while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Cell_line/){
        my $CELL_Line =$f[0];
        my $Histology =$f[1];
        $CELL_Line =~ s/\_//g;
        $CELL_Line =~ s/\-//g;
        $CELL_Line =~ s/\;//g;
        $CELL_Line =~ s/\,//g;
        $CELL_Line =~ s/\s+//g;
        $CELL_Line = uc($CELL_Line);
        push @{$hash2{$CELL_Line}},$Histology;
    }
}

while(<$I2>){
    chomp;
    my @f = split/\t/;
    my $CELL_Line = $f[0];
    my $disease = "Cancer";
    $CELL_Line =~ s/\_//g;
    $CELL_Line =~ s/\-//g;
    $CELL_Line =~ s/\;//g;
    $CELL_Line =~ s/\,//g;
    $CELL_Line =~ s/\s+//g;
    $CELL_Line = uc($CELL_Line);
    # push @{$hash2{$Cell_line}},$disease;
    push @{$hash2{$CELL_Line}},$disease;
}


while(<$I3>){
    chomp;
    my @f = split/\t/;
    unless(/^CELL_Line/){
        my $CELL_Line =$f[0];
        my $disease = $f[1];
        $CELL_Line =~ s/\_//g;
        $CELL_Line =~ s/\-//g;
        $CELL_Line =~ s/\;//g;
        $CELL_Line =~ s/\,//g;
        $CELL_Line =~ s/\s+//g;
        $CELL_Line = uc($CELL_Line);
        push @{$hash2{$CELL_Line}},$disease;
        # if($CELL_Line =~/GM12878/){
        #     print "$CELL_Line\n";
        # }
    }
}





while(<$I4>){
    chomp;
    my @f = split/\s+/;
    my $link = $f[2];
    my @t = split/\//,$link;
    # $link =~ s///g;
    # print  "$t[7]\n"; 
    my $file_name = $t[7];
    $file_name =~ s/wgEncodeUwTfbs//g;   
    $file_name =~ s/Ctcf.*//g;   
    # wgEncodeUwTfbsA549CtcfStdPkRep2.narrowPeak.gz
    unless(exists $hash1{$file_name}){
        $hash1{$file_name} =1;
        print $O1 "$file_name\n";
    }
    my $CELL_Line = $file_name;
    $CELL_Line =~ s/\_//g;
    $CELL_Line =~ s/\-//g;
    $CELL_Line =~ s/\;//g;
    $CELL_Line =~ s/\,//g;
    $CELL_Line =~ s/\s+//g;
    $CELL_Line = uc($CELL_Line);
    if (exists $hash2{$CELL_Line}){
        my @vs = @{$hash2{$CELL_Line}};
        my $disease = join("|",@vs);
        print $O2 "$file_name\t$disease\n";
    }
    else{
        print $O3 "$file_name\n";
    }
    
}




