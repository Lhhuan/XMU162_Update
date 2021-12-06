# for j in `seq 4044786 4044815`;do
# 	for i in `awk 'NR>1{print $0}' /home/liuke/project/single.cell/uterus/ATAC/markgenes.info.csv|cat`;do
# 		chr=$(echo $i |cut -f1 -d ",")
# 		star=$(echo $i |cut -f5 -d ",")
# 		end=$(echo $i |cut -f6 -d ",")
# 		gene=$(echo $i |cut -f4 -d ",")
# 		awk -v chr=$chr -v star=$star -v end=$end '$1==chr && $2>=star && $3<=end{print $0}' /home/liuke/project/single.cell/uterus/cutaneous_T_cell_lymphoma/SRR$j/SRR$j.bed > /home/liuke/project/single.cell/uterus/ATAC/csv/SRR$j$gene.csv
# 	done
# done


# for i in `awk -F"," 'NR>1{print $4}' /home/liuke/project/single.cell/uterus/ATAC/markgenes.info.csv|cat`;do
# 	cat /home/liuke/project/single.cell/uterus/ATAC/csv/*$i.csv > /home/liuke/project/single.cell/uterus/ATAC/bed/$i.bed
# done

# for i in `ls *.bed`;do
# 	sed -i 's/chr//g' $i
# done


# for i in `awk -F"," 'NR>1{print $7}' markgenes.info.csv |sort|uniq`;do
# 	cluster=$i
# 	for j in `awk -F"," -v cluster=$cluster '$7==cluster{print $4}' /home/liuke/project/single.cell/uterus/ATAC/markgenes.info.csv`;do
# 		touch /home/liuke/project/single.cell/uterus/ATAC/bed/$i.bed
# 		cat /home/liuke/project/single.cell/uterus/ATAC/bed/$j.bed >> /home/liuke/project/single.cell/uterus/ATAC/bed/$i.bed
# 	done
# done




for i in `ls cluster*.bed|tail -14`;do
	name=$(echo $i|cut -f1 -d ".")
	cat /home/liuke/project/single.cell/uterus/ATAC/joo > jobs
	echo "findMotifsGenome.pl /home/liuke/project/single.cell/uterus/ATAC/bed/$i /home/liuke/reference/Genome/hg19/hg19.fa /home/liuke/project/single.cell/uterus/ATAC/motif/${name}_motif -len 8,10,12 -size 200 -preparsedDir /home/liuke/reference/homer" >> jobs
	qsub jobs
done
