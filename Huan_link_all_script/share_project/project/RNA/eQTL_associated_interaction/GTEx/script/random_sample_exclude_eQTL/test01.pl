#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

use strict;
use warnings;
use Data::Dumper;

#  #!/usr/bin/perl
# my %hash;
# my @array = (1,2,3,2,4,5,5);

# foreach my $aaa(@array){
#     my $bb =  ++$hash{$aaa};
#     print "$bb\n";
# }



# # #grep 保存符合条件的元素
# # @array = grep { ++$hash{$_} >1 } @array;

# # #print join(" ",@array);
# # foreach my $item (@array){
# #     print $item;
# #     print "\n";
# # }


# #!/usr/bin/perl
# use strict;
# use warnings;
# use Data::Dumper;

# my @a = (1,2,3,4,5,6,7,8);
# my @b = (5,6,7,8,9,10);
# my %hash_a = map{$_=>1} @a;
# # my %hash_b = map{$_=>1} @b;
# # my %merge_all = map {$_ => 1} @a,@b;
# # my @a_only = grep {!$hash_b{$_}} @a;
# my @b_only = grep {!$hash_a{$_}} @b;
# # print Dumper(\@b_only);
# print "@b_only\n";




# my $dna='AACCGTTAATGGGCATCGATGCTATGCGAGCT';
# srand(time|$$);
# for (my $i=0;$i<20;++$i)
# {
# 	print randomposition($dna), " ";
# }
# print "\n";
# exit;

# sub randomposition
# {
# 	my($string)=@_;
# 	return int rand length $string;
# }

#--------------

# #!/usr/bin/perl
# use strict;
# use warnings;

# my $size=3;
# my $max=30;
# my $min=20;

# my @random=( );
# srand(time|$$);


# my $randna=dnaset($min,$max,$size);
# #print $randna;
# =pod
# foreach my $dna (@randna) {
#  print "$dna\n";
# }
# print "\n";
# exit;
# =cut

# #随机DNA的数据集。
# sub dnaset{
# my @set;
# for (my $i=0;$i<$size ;$i++) {
#     my($length)=randlength($min,$max);
# #获取随机长度。
# my($dna)=makedna($length);
# push(@set,$dna);
# print $dna,"\n";
# #用push在后面加入。
#     }
# return @set;
# }

# #用最大长度及最小长度，来控制，以获得一个随机长度。
# #返回的是一个整数值，最小值作为一个底数，
# #底数上加一个由最大值和最小值控制的随机数。
# sub randlength{
#     my($min,$max)=@_;
# return (int(rand($max-$min+1))+$min);
# }

# #把获得的随机碱基串起来。
# sub makedna{
# my($length)=@_;
# my $dna;
# for (my $i=0;$i<$length ;$i++) {
#          $dna.=randnucle( );
#     }
# return $dna;
# }

# #获取随机的碱基。
# sub randnucle{
# my (@nucles)=('A','C','G','T');
# #return  randele[rand @nucles];
#     return $nucles[rand @nucles];
# }




#-------------
my $dna='AACCGTTAATGGGCATCGATGCTATGCGAGCT';
srand(time|$$);
for (my $i=0;$i<20;++$i)
{
	print randomposition($dna), " ";
}
print "\n";
exit;

sub randomposition
{
	my($string)=@_;
	return int rand length $string;
}