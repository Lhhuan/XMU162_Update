perl 01_random_select_hotspot.pl
echo "01 random\n"
perl 07_merge_hic_left_prediction.pl
echo "07 left\n"
perl 07_merge_hic_right_prediction.pl
echo "07 right\n"
perl 08_merge_left_rigth_overlap.pl
echo "08 merge\n"