#把../output/07_logfold_merge_drug_ccle.txt和../data/secondary-screen-replicate-collapsed-treatment-info.txt
#merge 到一起，并把../data/secondary-screen-replicate-collapsed-treatment-info.txt中的drug 的disease area merge 到一起
#得../output/12_merge_drug_class_and_log_fold_change.txt.gz，并将log-fold change 取median得../output/12_merge_drug_class_and_log_fold_change_median.txt.gz
#
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use List::Util qw/max min/;

# my $f1 = "../data/secondary-screen-replicate-collapsed-treatment-info.txt";
my $f1 = "secondary-screen-dose-response-curve-parameters.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# my $f2 = "../output/07_logfold_merge_drug_ccle.txt.gz";
# open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
# open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
# my $fo1 = "../output/12_merge_drug_class_and_log_fold_change.txt.gz"; #此时drug cancer pair 的auc是unique的
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
# my $fo2 = "../output/12_merge_drug_class_and_log_fold_change_median.txt.gz";
# open my $O2, "| gzip >$fo2" or die $!;
# open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
# my $fo3 = "../output/13_prediction_and_Depmap_X3.txt";
# open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";

my (%hash1,%hash2,%hash3,%hash4,%hash5,%hash6,%hash7);




while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless (/^column_name/){
        my $column_name = $f[0]; #broad_id
        my $broad_id =$f[0];
        my $name =$f[12];
        my $dose ="NA";
        my $screen_id = $f[3];
        my $disease_area =$f[-6];
        $column_name =~ s/"//g;
        $disease_area =~ s/"//g;
        $dose =~ s/"//g;
        $broad_id =~ s/"//g;
        $name =~ s/"//g;
        $screen_id =~ s/"//g;
        my @ts = split/,/,$disease_area;
        foreach my $t(@ts){
            $t =~s/^\s+//;
            print "$t\n";
        }
    }
}


# while(<$I2>)
# {
#     chomp;
#     my @f= split/\t/;
#     if (/^broad_id/){
#         # print $O1 "$_\tdisease.area\tdrug_type\n";
#         print $O1 "$_\tdisease.area\tdrug_type\tbroad_id\tname\tdose\tscreen_id\n";
#     }
#     else{
#         my $Drug_chembl_id_Drug_claim_primary_name = $f[2];
#         my $drug =$f[0]; #broad_id
#         my $Sample = $f[4]; #DepMap_ID
#         my $Value =$f[1]; #log2_fold_change
#         my $k = "$Drug_chembl_id_Drug_claim_primary_name\t$Sample";
#         if(exists $hash1{$drug}){  #
#             unless($Value =~/NA/){
#                 my $drug_type = $hash1{$drug};
#                 print $O1 "$_\t$drug_type\n";
#                 push @{$hash2{$k}},$Value;
#             }
#         }
#     }
# }

# close($O1);

# #-------------
# my $f3 =$fo1;
# open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

# while(<$I3>)
# {
#     chomp;
#     my @f= split/\t/;
#     unless (/^broad_id/){
#         my $drug =$f[0]; #broad_id
#         my $Value =$f[1]; #primary_replicate_collapsed_logfold_change_value
#         my $Drug_chembl_id_Drug_claim_primary_name = $f[2];
#         my $Sample = $f[4]; #DepMap_ID
#         my $disease_area =$f[-6];
#         my $drug_type =$f[-5];
#         my $k = "$Drug_chembl_id_Drug_claim_primary_name\t$Sample\t$Value";
#         my $v = "$drug\t$disease_area\t$drug_type";
#         $hash3{$k}=$v;
#     }
# }
# print $O2 "Drug_chembl_id_Drug_claim_primary_name\tDepMap_ID\tAUC\tbroad_id\tdisease_area\tdrug_type\n";

# foreach my $k2(sort keys %hash2){
#      my @vs = @{$hash2{$k2}};
#      my $number = @vs;
#     # #  print "$number\n";
#     #  if ($number >1 ){
#     #     #  print "$k2\t$number\n";
#     #      foreach my $v(@vs){
#     #          print "$k2\t$number\t$v\n";
#     #      }
#     #  }
#      my $median_value_in_paper = mid( @vs); #value_in_paper #
#      #   #We label cell lines as sensitive to a treatment if the median collapsed fold-change is less than 0.3.
#      # 此value 是表示viability values，越大表明cell viability 能力越强.文中说<0.3认为是敏感
#     #  my $k3 = "$k2\t$vs[0]"; #此构造只为取drug info
#      my $k3 = "$k2\t$median_value_in_paper"; 
#      if (exists $hash3 {$k3}){
#          my $v3 = $hash3 {$k3};
#         #  print $O2 "$k3\t$v3\n";
#          print $O2 "$k2\t$median_value_in_paper\t$v3\n";
#      }
#      else{
#          print "$k3\n";
#      }
# }



# #---------------------------------------------------------------------------- some function
# #------------------------------------------ median
# # sub mid{
# #     my @list = sort{$a<=>$b} @_;
# #     my $count = @list;
# #     if( $count == 0 )
# #     {
# #         return undef;
# #     }   
# #     if(($count%2)==1){
# #         return $list[int(($count-1)/2)];
# #     }
# #     elsif(($count%2)==0){
# #         return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
# #     }
# # }  
# sub mid{ #三个取中间，两个取最小
#     my @list = sort{$a<=>$b} @_;
#     my $count = @list;
#     if( $count == 1 )
#     {
#         return $list[0];
#     }   
#     elsif($count==2){
#         return $list[0];
#     }
#     elsif($count==3){
#         return $list[1];
#     }
# }  
# close ($O1);