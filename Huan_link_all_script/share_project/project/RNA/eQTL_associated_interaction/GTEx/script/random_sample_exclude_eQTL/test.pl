#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;

# my @a = (1,2,3,4,5,6,7,8);
# my @b = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);



# my %hash_a = map{$_=>1} @a;
# my %hash_b = map{$_=>1} @b;

# my @b_only = grep {!$hash_a{$_}} @b;
# foreach my $b(@b_only){
#     print "$b\n";
# }
my @SNP_poss = ("A".."ZZ");
my $QTL_number = 6;

# my $pm = Parallel::ForkManager->new(10);
for(my $i=1; $i<11;$i++){
    # my $pid = $pm->start and next; #开始多线程
    my @new_poss = randomElem ( $QTL_number, @SNP_poss ) ;
    my %hash;
    @new_poss = grep { ++$hash{$_} < 2 } @new_poss;
    my $output= join("\t",@new_poss);
    print "$output\n";
    my $number = @new_poss;
    print "$number\n";

    # $pm->finish;  #多线程结束
}


sub randomElem { #随机取
    my ($want, @array) = @_ ;
    my (%seen, @ret);
     while ( @ret != $want ) {
     my $num = abs(int(rand(@array))); #@array 是指数组的长度，而$#array是指最后一个索引，由于rand的特殊性，如果用$#array会导致取不到最后一个值。
      if ( ! $seen{$num} ) { 
        # ++$seen{$num};
        # my $aaa  = $seen{$num};
        # print "$aaa\n";
        push @ret, $array[$num];
      }
     }
    return @ret;     
}






