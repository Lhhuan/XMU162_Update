2015-02-12
Marina Lizio <marina.lizio@riken.jp>
This folder contains the updated set of CAGE based enhancers for human and mouse phase 1 and 2 combined data.
Files

*human_permissive_enhancers_phase_1_and_2.bed.gz
*mouse_permissive_enhancers_phase_1_and_2.bed.gz

contain the coordinates of the enhancers. Coordinates are relative to hg19 and mm9 assemblies for human and mouse respectively.
Each predicted enhancer is described in BED12 format (columns 1-6 as chromosome, start coordinate, end coordinate, enhancer ID, score and strand) with two blocks denoting the merged regions of transcription initiation on the minus and plus strands. The thickStart and thickEnd columns denote the inferred mid position. The score column gives the maximum pooled expression of TCs used to construct each bidirectional loci.

Files

*human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt.gz
*mouse_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt.gz
*human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz
*mouse_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz

contain expression matrices representing enhancers expression in all samples.
Expression is provided in both raw counts and RLE-normalized TPM (tags per million) counts format.
Data provided by Robin Andersson <robin@binf.ku.dk> and Albin Sandelin <albin@binf.ku.dk>

#updates August 10 2015
Added Sample name to sample ID tables for human and mouse.

