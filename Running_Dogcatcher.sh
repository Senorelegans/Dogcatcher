#!/bin/bash

GTF=/gtf/Homo_sapiens.GRCh38.90.gtf
IN=/Dogcatcher_output
COMPARISON_PATH=${IN}/Treatment_vs_Control
BAMPATH=BAMS_FOLDER
BEDPATH=BEDS_FOLDER
CPUS=4
PADJ=0.05


mkdir -p ${IN}
mkdir -p ${COMPARISON_PATH}



python 1.0_Dogcatcher_flatten_gtf.py --annotation_file_with_path ${GTF}


###############= For Single File ###################
python 2.0_Dogcatcher.py \
--cpus $CPUS \
--BedGraph_input_min_strand ${BEDPATH}/Treatment_1_min.BedGraph \
--BedGraph_input_plu_strand ${BEDPATH}/Treatment_1_plu.BedGraph \
--output_prefix  ${IN}/ \
--annotation_file_with_path ${GTF} \
--window_size 100 \
--coverage_percentage 90


############### If running multiple files
###############Loop through POS and MIN Files ###################
declare -a MIN=($(ls -d ${BEDPATH}/*min.BedGraph | sort))
declare -a POS=($(ls -d ${BEDPATH}/*plu.BedGraph | sort))
L=${#MIN[@]}
X=0
while [  $X -le $L ]; do
    echo "**********"
    echo ${MIN[X]}
    echo ${POS[X]}
    python 2.0_Dogcatcher.py \
    --cpus $CPUS \
    --BedGraph_input_min_strand ${MIN[X]} \
    --BedGraph_input_plu_strand ${POS[X]} \
    --output_prefix  ${IN}/ \
    --annotation_file_with_path ${GTF} \
    --window_size 200 \
    --coverage_percentage 80
    let X=X+1
done
################################################################



python 2.5_Dogcatcher_filter.py \
--filter longest \
--input_prefix ${IN}/ \
--Dogcatcher_plu_strand_list \
Treatment_1_plu.BedGraph \
Treatment_2_plu.BedGraph \
Treatment_3_plu.BedGraph \
Control_1_plu.BedGraph \
Control_2_plu.BedGraph \
Control_3_plu.BedGraph \
--Dogcatcher_min_strand_list \
Treatment_1_min.BedGraph \
Treatment_2_min.BedGraph \
Treatment_3_min.BedGraph \
Control_1_min.BedGraph \
Control_2_min.BedGraph \
Control_3_min.BedGraph \
--output_prefix ${COMPARISON_PATH}/

python 3.0_Create_R_subread_DESeq2_script.py \
--annotation_file_with_path ${GTF} \
--control_BAM_list \
${BAMPATH}/Control_1.bam \
${BAMPATH}/Control_2.bam \
${BAMPATH}/Control_3.bam \
--treatment_BAM_list \
${BAMPATH}/Treatment_1.bam \
${BAMPATH}/Treatment_2.bam \
${BAMPATH}/Treatment_3.bam \
--input_R_template_file R_subread_DEseq2_TEMPLATE.R \
--input_prefix ${IN} \
--output_prefix ${IN}/initial_Rsubread_DESeq2 \
--cpus ${CPUS} \
--padj ${PADJ}

###Run the created R script
R CMD BATCH ${IN}/initial_Rsubread_DESeq2/Rsubread_DESeq2_initial.R ${IN}/initial_Rsubread_DESeq2/Rsubread_DESeq2_initial.R.out

##Run the second python script to generate DOG gtf with non-significant genes for proper normalization
python 4.0_Dogcatcher_Rsubread_DESeq2.py \
--annotation_file_with_path ${GTF} \
--input_prefix ${COMPARISON_PATH} \
--input_prefix_DESeq2 ${IN}/initial_Rsubread_DESeq2 \
--output_prefix ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes \
--padj ${PADJ}


##Run the R scripts from read-through gtf's combined with non-significant genes
R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig.R.out
R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig.R.out
R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsig.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsig.R.out
R CMD BATCH ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsig.R ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsig.R.out
echo finished R scripts

#Filter out significant read-through from DESeq2 and match with the csv so you can get biotypes etc.
python 5.0_filter_sig_DESeq2.py \
--annotation_file_with_path ${GTF} \
--input_prefix ${COMPARISON_PATH} \
--output_prefix ${COMPARISON_PATH}/FINAL_OUT \
--input_DOG_DESeq2_plu_sense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv \
--input_DOG_DESeq2_min_sense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_DOG_with_nonsig_Dogcatcher_DESeq2_sense.csv \
--input_ADOG_DESeq2_plu_antisense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/plu_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv \
--input_ADOG_DESeq2_min_antisense_file ${COMPARISON_PATH}/Dogcatcher_with_non-significant_genes/min_ALL_SAMPLES_ADOG_with_nonsig_Dogcatcher_DESeq2_antisense.csv \
--padj ${PADJ}
wait
