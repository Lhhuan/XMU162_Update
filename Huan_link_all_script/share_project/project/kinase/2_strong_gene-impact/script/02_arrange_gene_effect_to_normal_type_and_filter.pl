# 把../data/DepMap/Achilles_gene_effect.csv转换成常用的格式，即，gene cell_line probability
#得../output/02_arrange_gene_effect_to_normal_type.txt,
#过滤得dependency probability>0的文件得../output/02_arrange_gene_effect_to_normal_type_dependency.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../data/DepMap/Achilles_gene_effect.csv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "../output/02_arrange_gene_effect_to_normal_type.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/02_arrange_gene_effect_to_normal_type_dependency.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my (%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7,%hash8,%hash9,%hash10,%hash11,%hash12,%hash13,%hash30);

my $header = "cell_lines(Broad_IDs)\tHUGO(Entrez)\tCERES";
print $O1 "$header\n";
print $O2 "$header\n";

my @samples =();
while(<$I1>)
{
    chomp;
    my @f= split/\,/;
    if ($. ==1){
        push @samples,@f;#将sample 存进数组
        # my $i =0;
        # foreach my $title(@samples){
        #     $i++;
        #     print STDERR "$i\n"; #可以看出，一共5449列，
        # }
        # print STDERR "$samples[5448]\n";
    }
    else{
        for (my $i=1;$i<18334;$i++){ #对原表格进行逐个转换
            my $cell_line = $f[0];
            my $gene = $samples[$i];
            my $probability = $f[$i];
            my $output = "$cell_line\t$gene\t$probability";
            # print "$probability\n";
            print $O1 "$output\n";
            unless($probability =~/NA/){ #过滤dependency probability>0
                if ($probability < -0.6){ #
                    print $O2 "$output\n";
                }               
            }
        }  
    }
}

close($O1);