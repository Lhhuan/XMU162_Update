perl 01_merge_enhancer.pl #将./enhancer 下的文件合并起来得 01_merge_enahcner_sample.bed.gz 
liftOver 01_merge_enahcner_sample.bed.gz  /home/huanhuan/reference/hg19ToHg38.over.chain.gz ./hg38/01_merge_enhancer_sample.bed  01_unmap_merge_enahcner_sample.bed 
less ./hg38/01_merge_enhancer_sample.bed | sort -k1,1 -k2,2n |gzip >./hg38/01_merge_enhancer_sample_sorted.bed.gz
perl 02_merge_enhancer_target.pl #将./enhancer_gene 下的文件合并起来得 02_merge_enhancer_target_sample.bed.gz, sorted  02_merge_enhancer_target_sample_sorted.bed.gz
# python ~/.conda/envs/huan_py3/bin/CrossMap.py bed /home/huanhuan/reference/hg19ToHg38.over.chain.gz 02_merge_enhancer_target_sample.bed.gz ./hg38/02_merge_enhancer_target_sample.bed
# less ./hg38/02_merge_enhancer_target_sample.bed |sort -k1,1 -k2,2n |gzip >./hg38/02_merge_enhancer_target_sample.bed.gz