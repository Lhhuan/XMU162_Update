##normalize "../raw_data/interactions.tsv.gz", 得../normalized/01_interactions_normalized.txt;
less ../normalized/01_interactions_normalized.txt |cut -f6,15 |sort -u | wc -l #1806143 -1