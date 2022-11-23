#选出normal cell 的POL2

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $f1 = "../cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 



my %hash1;
while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $CCLE_name = $f[0];
    my $Disease = $f[1];
    unless(/^CELL_Line/){
        if ($Disease =~/NO/){
            $hash1{$CCLE_name}=1;
            # print "$CCLE_name\n";
        }
        
    }
}


my %hash5;
for (my $i=1;$i<23;$i++){
    my $k = "chr${i}";
    $hash5{$k}=1;
}


# foreach my $dir(@dirs){
    my (%hash2,%hash3,%hash4);
    #----------------read QC file;
    my $dir ="../Human_FACTOR/human_factor";
    my $f2 = "${dir}_full_QC.txt";
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  

    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        my $DCid = $f[0];
        my $GSMID = $f[2];
        my $factor = $f[3];
        my $Cell_line = $f[4];
        my $v = "$factor\t$Cell_line";
        if (exists $hash1{$Cell_line}){
            if($factor =~/\bPOLR2A|POLR2AphosphoS5|POLR2B|POLR2H|POLR2M\b/){
                $hash2{$DCid}=$v;
                # print "$factor\n";
            }
            
        }
    }

#--------------------------------
    opendir (DIR, $dir) or die "can't open the directory!";
    my @files = readdir DIR;
    #---------------------------------output option
    my @aa = split/\//,$dir;
    my $out_dir = $aa[1];
    # print "$out_dir\n";
    my $output_dir = "../normal_cell/${out_dir}/POL2";
    my $fo1 = "${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "${output_dir}/merge_pos_info_narrow_peak.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    my $fo3 = "${output_dir}/merge_pos_info_narrow_peak_signalValue.bed.gz";
    open my $O3, "| gzip >$fo3" or die $!;
    # print $O1 "chr\tstart\tend\tfactor\tcell_line\tfile_name\n";
    #--------对dir 下all file 进行判断
    foreach my $file(@files){
        if ( $file =~ /[a-z]/) {
            my @t = split/\_/,$file;
            my $id = $t[0];
            # print "$id\n";
            if (exists $hash2{$id}){
                my $v =$hash2{$id};
                my @h = split/\t/,$v;
                my $factor = $h[0];
                my $cell_line =$h[1];
                # print "$id\n";
                my $f3 = "$dir/$file";
                # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
                open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
                # print "$f3\n";
                while(<$I3>)
                {
                    chomp;
                    my @f = split/\s+/;
                    my $chr = $f[0];
                    my $start =$f[1];
                    my $end =$f[2];
                    my $name = $f[3];
                    my $signalValue =$f[6];
                    if (exists $hash5{$chr}){
                        my $output1 = "$chr\t$start\t$end\t$factor\t$cell_line\t$file";
                        my $output2 = "$chr\t$start\t$end";
                        my $output3 = "$chr\t$start\t$end\t$signalValue";
                        print $O3 "$output3\n";
                        unless(exists $hash3{$output1}){
                            $hash3{$output1}=1;
                            print $O1 "$output1\n";
                        }
                        unless(exists $hash4{$output2}){
                            $hash4{$output2}=1;
                            print $O2 "$output2\n";
                        }
                    }
                }
            }
        }
    }
    close($O1);
    close($O3);
    system "zless $fo1 |sort -k1,1 -k2,2n |gzip > ${output_dir}/merge_pos_info_sample_narrow_peak_sorted.bed.gz";
    system "zless $fo2 |sort -k1,1 -k2,2n |gzip > ${output_dir}/merge_pos_info_narrow_peak_sorted.bed.gz"
    system "zless $fo3 |sort -k1,1 -k2,2n |gzip > ${output_dir}/merge_pos_info_narrow_peak_signalValue_sorted.bed.gz";
# }
