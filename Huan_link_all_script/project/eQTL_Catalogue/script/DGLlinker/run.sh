Rscript adjust_Inact_format.R #adjust_gene_IntAct.csv
echo "Source node name,Source node type,Relationship type,Target node type,Target node name" > ./train/new_adjust_gene_IntAct.csv
sed -n '1!p' ./train/adjust_gene_IntAct.csv >> ./train/new_adjust_gene_IntAct.csv
python train_validation_models.py