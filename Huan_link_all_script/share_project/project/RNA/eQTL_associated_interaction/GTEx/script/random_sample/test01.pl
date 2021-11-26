#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

use strict;
use warnings;
use Data::Dumper;

 #!/usr/bin/perl
my %hash;
my @array = (1,2,3,2,4,5,5);

foreach my $aaa(@array){
    my $bb =  ++$hash{$aaa};
    # print "$bb\n";
}

my $count= keys %hash;
print "$count\n";

# #grep 保存符合条件的元素
# @array = grep { ++$hash{$_} >1 } @array;

# #print join(" ",@array);
# foreach my $item (@array){
#     print $item;
#     print "\n";
# }