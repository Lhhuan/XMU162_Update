
cluster_inputdir="/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200"
bw_dir="/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/"
level="mean"
computeMatrix scale-regions \
            -S ${bw_dir}/H3K27ac_merge_${level}_signalvalue.bw  \
            -R ${cluster_inputdir}/cluster_1.bed \
                ${cluster_inputdir}/cluster_2.bed \
                ${cluster_inputdir}/cluster_3.bed \
                ${cluster_inputdir}/cluster_4.bed \
                ${cluster_inputdir}/cluster_5.bed \
                ${cluster_inputdir}/cluster_6.bed \
                ${cluster_inputdir}/cluster_7.bed \
            -out /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/test_${level}.tab.gz
            # --startLabel border\
            # --endLabel border
            # --skipZeros \

plotProfile -m /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/test_${level}.tab.gz \
                -out /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/test_region_${level}.pdf\
                --startLabel border\
                --endLabel border


plotHeatmap -m /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/test_${level}.tab.gz \
                -out /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/test_region_${level}_heatmap.pdf\
                --startLabel border\
                --endLabel border\
                # --colorMap "#42B540FF" "#0099B4FF" "#A593E0" "#fcbe32" "#F16B6F" "#F68657" "#5A9367"



plotProfile -m ${output_dir}/${marker}_${level}.tab.gz \
              -out ${output_dir}/${marker}_${level}.pdf\
              --colors "#42B540FF" "#0099B4FF" "#A593E0" "#fcbe32" "#F16B6F" "#F68657" "#5A9367"\
              --plotHeight 10.5 \
              --plotWidth 10 \
              --refPointLabel 'center'  --regionsLabel "C1" "C2" "C3" "C4" "C5" "C6" "C7"\
              --yAxisLabel "" 

