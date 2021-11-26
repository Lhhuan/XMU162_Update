use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;

my (%hash1,%hash2,%hash3,%hash4);
my @qtl=("cis_1MB","trans_1MB","cis_10MB","trans_10MB");
$hash1{kkk}=1;
$hash2{kkk}=2;
$hash3{kkk}=3;
$hash4{kkk}=4;
my @h = ({%hash1},{%hash2},{%hash3},{%hash4});

for (my $i=0;$i<4;$i++){
# print "$h[$i]{kkk}\n";
print "$h[$i]\n";
}