#We performed eQTL mapping in normal skin from healthy controls (normal skin, n=57), 
#normal-appearing skin from psoriatic patients (uninvolved skin, n=53), and diseased skin from psoriatic patients (lesional skin, n=53).
wget -c http://csg.sph.umich.edu/junding/eQTL/TableDownload/NN57subj_p1E-5_annotate_allcis_1Mb.tbl
mv NN57subj_p1E-5_annotate_allcis_1Mb.tbl significant_cis-eQTLs_from_normal_skin.txt
#---------------
wget -c http://csg.sph.umich.edu/junding/eQTL/TableDownload/PN53subj_p1E-5_annotate_allcis_1Mb.tbl
mv PN53subj_p1E-5_annotate_allcis_1Mb.tbl significant_cis-eQTLs_from_uninvolved_skin.txt
#----------
wget -c http://csg.sph.umich.edu/junding/eQTL/TableDownload/PP53subj_p1E-5_annotate_allcis_1Mb.tbl
mv PP53subj_p1E-5_annotate_allcis_1Mb.tbl significant_cis-eQTLs_from_lesional_skin.txt