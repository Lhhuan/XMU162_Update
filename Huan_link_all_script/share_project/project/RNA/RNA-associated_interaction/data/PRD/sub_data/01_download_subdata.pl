#download data 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

for (my $i=1;$i<1083;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
my $link = "http://pri.hgc.jp/export?p=${i}&c=10";
my $commod = "wget ${link}\n";
system "$commod";
}