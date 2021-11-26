
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use JSON;
use Encode;
use Data::Dumper;

# my (%hash1, %hash2, %hash3,%hash5);
# my $f1 = "../data/QTLbase_download_data_sourceid.txt";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "sub_table.json";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
# my $fo1 = "../output/chr22.txt";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# my $fo2 = "../output/chr22_eur.txt";
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

my ($json) = @ARGV;

my $context;

while (<$I2>) {
	$context .= $_;
    # chomp;
    # print "$_\n";
}
close $I2;

my $obj    = decode_json($context);
printf Dumper($obj)."\n";
# print "$obj\n";






# while(<$I1>)
# {
#     chomp;
#     unless(/PMID/){
#         # print "$file\n";
#         my @f =split/\t/;
#         my $Sourceid = $f[1];
#         my $Tissue =$f[3];
#         my $Population = $f[4];
#         my $v = "$Tissue\t$Population";
#         $hash5{$Sourceid}=$v; #用于为$f2添加$Tissue和$Population
#     }
# }


# while(<$I2>)
# {
#     chomp;
#     if(/SNP_chr/){
#         print $O1 "$_\tTissue\tPopulation\n";
#         print $O2 "$_\tTissue\tPopulation\n";
#     }
#     else{
#         my @f =split/\t/;
#         my $snp_chr = $f[0];
#         my $snp_pos = $f[1];
#         my $Sourceid =$f[-2];
#         if ($snp_chr =~/\b22\b/){
#             if (exists $hash5{$Sourceid}){
#                 my $v = $hash5{$Sourceid};
#                 print $O1 "$_\t$v\n";
#                 my @t =split/\t/,$v;
#                 my $pop = $t[1];
#                 if ($pop =~/EUR/){
#                     print $O2 "$_\t$v\n";
#                 }
#             }
#         }
#     }
# }