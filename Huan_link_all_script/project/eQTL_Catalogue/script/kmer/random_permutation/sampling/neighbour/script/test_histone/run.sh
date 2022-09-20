zless "/share/data0/GTEx/annotation/ROADMAP/sample/E019/E019-H3K27ac.narrowPeak.gz" >test.bed
# zless "/share/data0/GTEx/annotation/ROADMAP/sample/E020/E020-H3K27ac.narrowPeak.gz" >>test.bed
less test.bed| sort -k1,1 -k2,2n |cut -f1-3,7 >test_sort.bedgraph
cut -f1,2 "/share/Projects/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa.fai" >GRCh37_chrom.sizes
source activate huan_py3
bedGraphToBigWig test_sort.bedgraph GRCh37_chrom.sizes myBigWig.bw

conda deactivate


computeMatrix reference-point \
       -S myBigWig.bw  \
       -R /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200/cluster_3.bed \
          /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200/cluster_5.bed \
       --samplesLabel "NRF2 (A549)"\
       --referencePoint center \
       --binSize 50\
       -b 3000 -a 3000 \
       --skipZeros \
       -out nrf2_hepg2_sig1_3kb.tab.gz \


plotProfile -m nrf2_hepg2_sig1_3kb.tab.gz \
              -out nrf2_hepg2_sig1_3kb.pdf\
              --colors "#42B540FF" "#0099B4FF" \
              --plotHeight 4 \
              --plotWidth 4 \
              --refPointLabel 'center'  --regionsLabel "TSS" "Body"\
              --yAxisLabel "" \
            #   --xAxisLabel "Hypo MS1" \
            #   --plotTitle ""



#===========CTCF
/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38/normal_cell_line_ctcf_merge_mean_signalvalue.bw
normal_cell_line_ctcf_merge_median_signalvalue.bw
normal_cell_line_ctcf_merge_max_signalvalue.bw

#=======================cistrome
markers = c("Human_FACTOR","HISTONE_MARK_AND_VARIANT","Human_CHROMATIN_Accessibility")

/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/
marker/merge_mean_signalvalue.bw
marker/merge_median_signalvalue.bw
marker/merge_max_signalvalue.bw



#=====================histone
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3")
/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/
marker _merge_mean_signalvalue.bw
marker _merge_median_signalvalue.bw
marker _merge_max_signalvalue.bw