#将chr${i}.phastCons100way.wigFix.gz normalized成正常的文件，得./normalized_per_chr/chr${i}.normalized_merge_phastCons100way.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


# my $fo1 = "./merge_phastCons100way.bed.gz";
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;

# print $O1 "chr\tstart\tend\tvalue\n";
for (my $i=1;$i<23;$i++){
    my $input_name = "chr${i}.phastCons100way.wigFix.gz";
    my $f1 = $input_name;
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $fo1 = "./normalized_per_chr/chr${i}.normalized_merge_phastCons100way.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    my @starts=();
    my @steps = ();
    my @values=();
    while(<$I1>)
    {
        chomp;
        if (/^fixedStep/){
            @starts=();
            @steps = ();
            @values=();
            my @t = split/\s+/;
            my $start = $t[2];
            my $step =$t[3];
            $start =~ s/start=//g;
            $step =~ s/step=//g;
            # print "$step\t$i\n";
            push @starts, $start;
            push @steps,$step;
        }
        else{
            my $number = @values;
            if ($number <1){ #---------第一个是start
                my $start = $starts[0];
                my $end= $start +1;
                print $O1 "chr$i\t$start\t$end\t$_\n"; #the original data from ucsc, 0-bsed data, so the convert to bed file (start = end -1)
                push @values,$_;
            }
            else{ #-------------从第二个往后，$pos= $start +1
                my $step = $steps[0];
                my $start = $starts[0];
                my $pos = $start +$step;
                my $new_start = $pos;
                my $end = $new_start +1;
                print $O1 "chr${i}\t$new_start\t$end\t$_\n";
                @starts=();
                push @starts, $pos; #为一下计算做准备
            }
        }
    }
}