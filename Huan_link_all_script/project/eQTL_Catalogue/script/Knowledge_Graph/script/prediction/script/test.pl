#为./output/01_hotspot_target_gene_reactomeFI.bed.gz  annotation co-expression数据，得./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}，得./output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

my @f= ();
my $t=join(";",@f);
print "$t\n";