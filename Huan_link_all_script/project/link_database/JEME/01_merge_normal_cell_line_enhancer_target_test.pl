#将./enhancer 下的文件合并起来得 01_merge_enahcner_sample.bed.gz 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my (%hash1,%hash2);

my $f1 = "ENCODE_Roadmap_info_man.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $fo1 = "01_merge_enhancer_target_sample.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
while(<$I1>)
{
    chomp;
    unless(/^ID/){
        my @f = split/\t/;
        for (my $i=0;$i<8;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
            unless(defined $f[$i]){$f[$i] = "NA"};
        }
        my $ID = $f[0];
        my $Sample_name =$f[1];
        my $cancer =$f[-1];
        unless($cancer =~/YES/){
            $hash1{$ID}=$Sample_name;
        }
    }
}

my $f2 = "famtom5_info_man.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

while(<$I2>)
{
    chomp;
    unless(/^ID/){
        my @f = split/\t/;
        for (my $i=0;$i<7;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
            unless(defined $f[$i]){$f[$i] = "NA"};
        }
        my $ID = $f[0];
        my $Sample_name =$f[1];
        my $cancer =$f[-1];
        unless($cancer =~/YES/){
            $hash1{$ID}=$Sample_name;
        }
    }
}

my @dirs= ("fantom5_elasticnet","fantom5_lasso","encoderoadmap_lasso","encoderoadmap_elasticnet");
foreach my $dir(@dirs){
    print "$dir\tstart\n";
    opendir (DIR, $dir) or die "can't open the directory!";
    my @files = readdir DIR;
    foreach my $file(@files){
        if ( $file =~ /csv/) {
            my @t=split/\./,$file;
            my $ID =$t[1];
            if(exists $hash1{$ID}){
                my $Sample_name= $hash1{$ID};
                my $f3 = "$dir/$file";
                open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
                # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
                # print "$f3\n";
                while(<$I3>)
                {
                    chomp;
                    my @f = split/\,/;
                    my $pos= $f[0];
                    my $gene= $f[1];
                    my $score =$f[2];
                    $pos =~ s/:|-/\t/g;
                    # print "$f[0]\t$pos\n";
                    my @ff = split/\$/,$gene;
                    my $ensg=$ff[0];
                    if($ensg=~/ENSG/){
                        my $gene_symbol =$ff[1];
                        $ensg=~ s/\..*//g;
                        # print "$ff[0]\t$ensg\n";
                        print $O1 "$pos\t$ensg\t$Sample_name\t$dir\n";
                    }
                }               
            }
        }
    }
    print "$dir\tend\n";
}

#---------------------------------output option

#--------对dir 下all file 进行判断

close($O1);
system "zless $fo1 |sort -k1,1 -k2,2n |gzip > 01_merge_enhancer_target_sample_sorted.bed.gz";
