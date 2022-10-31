
cluster_inputdir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200"
bw_dir="/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/"
level="max"
marker="H3K36me3"
i=6
output_dir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/homer_200/test"
computeMatrix reference-point \
            -S ${bw_dir}/${marker}_merge_${level}_signalvalue.bw  \
            -R /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/homer_200/cluster_${i}.bed \
            --samplesLabel ${marker}\
            --binSize 50\
            -b 2500 -a 2500 \
            --referencePoint center \
            -out ${output_dir}/${i}_${marker}_${level}.gz

plotProfile -m ${output_dir}/${i}_${marker}_${level}.gz \
              -out ${output_dir}/${i}_${marker}_${level}.pdf\
              --plotHeight 4.5 \
              --plotWidth 10 \
              --refPointLabel 'center'  --regionsLabel $i\
              --yAxisLabel "" 



#------------------------------------------------         

cluster_inputdir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/re_2e_5/homer_200/"
bw_dir="/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/"
level="max"
marker="H3K36me3"
i=8
output_dir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/re_2e_5/homer_200/test/"
computeMatrix reference-point \
            -S ${bw_dir}/${marker}_merge_${level}_signalvalue.bw  \
            -R ${cluster_inputdir}/cluster_${i}.bed \
            --samplesLabel ${marker}\
            --binSize 50\
            -b 2500 -a 2500 \
            --referencePoint center \
            -out ${output_dir}/${i}_${marker}_${level}.gz

plotProfile -m ${output_dir}/${i}_${marker}_${level}.gz \
              -out ${output_dir}/${i}_${marker}_${level}.pdf\
              --plotHeight 4.5 \
              --plotWidth 10 \
              --refPointLabel 'center'  --regionsLabel $i\
              --yAxisLabel "" 
