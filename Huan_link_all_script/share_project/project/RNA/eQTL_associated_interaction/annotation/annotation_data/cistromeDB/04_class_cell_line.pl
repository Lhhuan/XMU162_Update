#筛选@files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt") 中在"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CCLE/CCLE_sample_info_file_2012-10-18.txt"的文件，得存在的./cell_line_info/04_existing_ccle_cell_line_info.txt,得不存在但是cancer cell line ./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt，得附加原始信息文件 ./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique_info.txt
#得不在ccle中存在，从名字中无法判断cancer cell line的 unique cell line 得./cell_line_info/04_unique_cell_line_without_info.txt 得附加原始信息文件 ./cell_line_info/04_cell_line_without_info_original.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;



my @files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt");
# my $pm = Parallel::ForkManager->new(4); ## 设置最大的线程数目

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CCLE/CCLE_sample_info_file_2012-10-18.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "./cell_line_info/04_existing_ccle_cell_line_info.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;

my $fo4 = "./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt";
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";
# open my $O2, "| gzip >$fo2" or die $!;

my $fo5 = "./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique_info.txt";
open my $O5, '>', $fo5 or die "$0 : failed to open output file '$fo5' : $!\n";
# open my $O2, "| gzip >$fo2" or die $!;

my $fo2 = "./cell_line_info/04_unique_cell_line_without_info.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
# open my $O2, "| gzip >$fo2" or die $!;

my $fo3 = "./cell_line_info/04_cell_line_without_info_original.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
# open my $O2, "| gzip >$fo2" or die $!;


print $O1 "Cell_line\tSite_primary\tHistology\tHist_Subtype1\n";
print $O3 "Cell_line\tCell_type\tTissue_type\n";
print $O5 "Cell_line\tCell_type\tTissue_type\n";


my (%hash1,%hash2,%hash3,%hash4,%hash5);
while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $CCLE_name = $f[0];
    my $Site_primary = $f[4];
    my $Histology = $f[5];
    my $Hist_Subtype1 =$f[6];
    my $v= "$Site_primary\t$Histology\t$Hist_Subtype1";
    unless(/^CCLE/){
        my @t = split/_/,$CCLE_name;
        my $cell_name_r = $t[0];
        $cell_name_r =uc($cell_name_r);
        $hash1{$cell_name_r}=$v;
    }
}


foreach my $file(@files){
    my $f2 = $file;
    # my $f2 = "./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
    # open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件



    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        my $DCid = $f[0];
        my $GSMID = $f[2];
        my $Cell_line = $f[4];
        my $Cell_type = $f[5];
        my $Tissue_type = $f[6];
        unless(/^DCid/){
            my $Cell_line_r = $Cell_line;
            $Cell_line_r =~ s/\s+//g;
            $Cell_line_r =~ s/\-//g;
            $Cell_line_r =~ s/\_//g;
            $Cell_line_r =uc($Cell_line_r);
            if (exists $hash1{$Cell_line_r}){
                my $v = $hash1{$Cell_line_r};
                # print $O1 "$Cell_line\t$Histology\n";
                my $output = "$Cell_line\t$v";
                unless(exists $hash3{$output}){
                    $hash3{$output}=1;
                    print $O1 "$output\n";
                }
            }
            else{
                unless($Cell_type =~ /\bNone\b/){
                    if ($Cell_type =~ /cancer|oma|Patient|tumor/i){
                        unless(exists $hash4{$Cell_line}){
                            $hash4{$Cell_line}=1;
                            print $O4 "$Cell_line\n";
                        } 
                        #---------------------
                        my $output5 = "$Cell_line\t$Cell_type\t$Tissue_type";
                        unless(exists $hash5{$output5}){
                            $hash5{$output5}=1;
                            print $O5 "$output5\n";
                        }                  
                    }
                    else{
                        unless(exists $hash2{$Cell_line}){
                            $hash2{$Cell_line}=1;
                            print $O2 "$Cell_line\n";
                        }
                        my $output2 = "$Cell_line\t$Cell_type\t$Tissue_type";
                        unless(exists $hash3{$output2}){
                            $hash3{$output2}=1;
                            print $O3 "$output2\n";
                        }
                    }
                }

                # if($Cell_type =~ /\bNone\b/){
                #     print "$Cell_type\n";
                # }
            }
        }
    }
}

    # $pm->finish;  #多线程结束

# }