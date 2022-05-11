
perl 001_adjust_hotspot.pl
python 01_graph_gen_huan.py
python 02_random_walk.py ../output/01_train_graph.dgl 50
python 02_random_walk.py ../output/01_test_graph.dgl 50
python 02_random_walk.py ../output/01_val_graph.dgl 50