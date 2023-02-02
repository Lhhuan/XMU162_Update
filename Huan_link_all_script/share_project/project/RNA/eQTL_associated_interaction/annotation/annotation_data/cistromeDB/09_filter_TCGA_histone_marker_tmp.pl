#筛选@files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt") 中在"./cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt"的cell lien，并在相应文件夹提取出文件得${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz，得对于文件及peak信息文件"${output_dir}/merge_pos_info_narrow_peak.bed.gz"

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my (%hashaaa,%hash1);

# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
# foreach my $cancer(@cancers){
#     $hashaaa{$cancer}=1;
# }
my $f1 = "./cancer.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    for (my $i=0;$i<5;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
        unless(defined $f[$i]){
            $f[$i] = "NONE";
        }
        unless($f[$i]=~/\w/){$f[$i]="NONE"}
    }
    my $CCLE_name = $f[1];
    my $TCGA = $f[4];
    $TCGA =~ s/"//g;
    # if (exists $hashaaa{$TCGA}){
    unless($TCGA =~ /NONE|Map_to_TCGA_Abbreviation/){
        $hash1{$CCLE_name}=$TCGA;
    }
}

my %hash5;
for (my $i=1;$i<23;$i++){
    my $k = "chr${i}";
    $hash5{$k}=1;
}

my @need_markers=("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3");
my %hash6;
for my $need_marker (@need_markers){
    $hash6{$need_marker}=1;
}

my $dir = "./HISTONE_MARK_AND_VARIANT/human_hm";
my (%hash2,%hash3,%hash4);
#----------------read QC file;
my $f2 = "${dir}_full_QC.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  

while(<$I2>)
{
    chomp;
    my @f = split/\t/;
    my $DCid = $f[0];
    my $GSMID = $f[2];
    my $factor = $f[3];
    # if($factor =~/^H3K|^H4K|^H3ac|^H4ac/){
    if($factor =~/H3K27ac|H3K9ac|H3K36me3|H3K4me1|H3K4me3|H3K27me3|H3K9me3/){
        # "H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3"
        # print "$factor\n";
        $factor =~s/\s+//g;
        my @ts = split/,/,$factor;
        foreach my $t(@ts){
            if (exists $hash6{$t}){
                my $Cell_line = $f[4];
                my $v = "$factor\t$Cell_line";
                if (exists $hash1{$Cell_line}){
                    my $TCGA = $hash1{$Cell_line};
                    my $v3 = "$DCid\t$Cell_line";
                    my $k = "$TCGA\t$t";
                    push @{$hash2{$k}},$v3;
                }
            }
        }
    }
}
my @aa = split/\//,$dir;
my $out_dir = $aa[1];

my $fo1 = "./cancer_cell/histone_TCGA_mark1.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my %hash_exist;
print $O1 "TCGA\tfactor\tcell_line\tDcid\n";
foreach my $k(sort keys %hash2){
    my @tt =split/\t/,$k;
    my $TCGA= $tt[0];
    my $factor = $tt[1];

    my @vs= @{$hash2{$k}};
    my $final_out_dir= "./cancer_cell/${out_dir}/${TCGA}/${factor}";
    unless(-e $final_out_dir){
        system "mkdir -p $final_out_dir";
    }

    # # open my $O1, "| gzip >$fo1" or die $!;
    # my $fo2 = "${final_out_dir}/merge_pos_info_narrow_peak.bed";
    # # open my $O2, "| gzip >$fo2" or die $!;
    #  open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
    # my $fo3 = "${final_out_dir}/merge_pos_info_narrow_peak_signalvalue.bed.gz";
    # open my $O3, "| gzip >$fo3" or die $!;
    # # open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
    # # print $O1 "chr\tstart\tend\tfactor\tname\tfile_name\n";
    foreach my $v(@vs){
        my @t = split/\t/,$v;
        my $DCid = $t[0];
        # my $factor = $t[1];
        my $Cell_line =$t[1];
        my $file = "${DCid}_sort_peaks.narrowPeak.bed.gz";
        my $f3 = "$dir/$file";
        if(-e $f3){
            my $out_infos= "$k\t$Cell_line\t$DCid";
            unless(exists $hash_exist{$out_infos}){
                $hash_exist{$out_infos}=1;
                print "$out_infos\n";
                print $O1 "$out_infos\n";
            }

            # # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
            # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
            # while(<$I3>)
            # {
            #     chomp;
            #     my @f = split/\s+/;
            #     my $chr = $f[0];
            #     my $start =$f[1];
            #     my $end =$f[2];
            #     my $name = $f[3];
            #     my $signalValue =$f[6];
            #     if (exists $hash5{$chr}){
            #         my $output1 = "$chr\t$start\t$end\t$factor\t$name\t$file";
            #         my $output2 = "$chr\t$start\t$end";
            #         my $output3 = "$chr\t$start\t$end\t$signalValue";
            #         print $O3 "$output3\n";
            #         # unless(exists $hash3{$output1}){
            #         #     $hash3{$output1}=1;
            #         #     print $O1 "$output1\n";
            #         # }
            #         unless(exists $hash4{$output2}){
            #             $hash4{$output2}=1;
            #             print $O2 "$output2\n";
            #         }
            #     }
            # }
        }
    }
    # close($O2);
    # close($O3);
    print "print_mark_file_finish\n";
    #----------------
    # unless(-z $fo2){  #判断文件是不为空
    #     my $out1 = "${final_out_dir}/merge_pos_info_narrow_peak_sort.bed.gz";
    #     my $out2 = "${final_out_dir}/merge_pos_info_narrow_peak_sort_union.bed.gz";
    #     # my $out3 = "${final_out_dir}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";
    #     my $out3 = "${final_out_dir}/merge_pos_info_narrow_peak_signalvalue_sort.bed.gz";
    #     my $out4 = "${final_out_dir}/merge_pos_info_narrow_peak_sort_union_signalvalue.bed.gz";
    #     my $command1 = "less $fo2 |sort -k1,1 -k2,2n |gzip >$out1 ";
    #     my $command2 = "bedtools merge -i $out1 |gzip > $out2";
    #     my $command3 = "zless $fo3 |sort -k1,1 -k2,2n |gzip >$out3";
    #     my $command4 = "bedtools intersect -a $out2 -b $out3  -wa -wb |gzip >$out4";
    #     system "$command1";
    #     system "$command2";
    #     system "$command3";
    #     system "$command4";
    #     system "gzip $fo2";
    # }

    # if(-z $fo2){
    #     system "rm -r $final_out_dir";
    # }
}

