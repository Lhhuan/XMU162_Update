#将../output/clump/*.clumped cut 出snp和与其相关的SNP，得../output/final_clump/*.clumped
#得全部的final_clump文件 ../output/07_all_clump_result.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "../output/07_all_clump_result.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print  $O1 "SNP\tSNP2\tSource\n";

#-------------------------------------------------------------#获取组织名称
my $dir = "../output/clump";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = ".clumped"; #文件名后缀.clumped
foreach my $file(@files){
    # my $suffix = ".clumped";
    if ($file =~/$suffix/){ 
        # $file =~ s/$suffix//g; #提取前缀disease
        # my $tissue =$file;
        # print "$tissue\n";
        # my $command  = "cat $file | cut -f >../output/final_clump/$file"
        my $f1 = "../output/clump/$file";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
        my $fo2 = "../output/final_clump/$file";
        open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
        # print "$f1\n";
        while(<$I1>)
        {
            chomp;
            unless(/CHR/){
                # print "$file\n";
                my @f =split/\s+/;
                for (my $i=0;$i<13;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
                    unless(defined $f[$i]){
                    $f[$i] = "NA";
                    }
                unless($f[$i]=~/\w/){$f[$i]="NULL"} #对文件进行处理，把所有定义的没有字符的都替换成NULL
                }
                # print "$f[-1]\n";
                my $SNP = $f[3];
                my $SNP2 = $f[-1];
                $file =~ s/$suffix//g; #提取前缀disease
                print $O1 "$SNP\t$SNP2\t$file\n";
                print $O2 "$SNP\t$SNP2\n";
            }
        }
    }
}