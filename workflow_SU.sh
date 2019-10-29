#!bin/bash

start1=$SECONDS

FASTQ1=<FASTQ1 FILE>
FASTQ2=<FASTQ2 FILE>
REF_GENOME=<REFERENCE GENOME>

OUTPUT_DIR=<PATH TO OUTPUT FOLDER>
TMPDIR=<PATH TO TEMPORARY DIRECTORY>

FGBIO=<PATH TO FGBIO>
PICARD=<PATH TO PICARD>
GATK=<PATH TO GATK>
BCFTOOLS=<PATH TO BCFTOOLS>
SCRIPTS=<PATH TO PYTHON SCRIPTS>

SMPL=<SAMPLE NAME>
NUCL=<HOW MANY NUCLEOTIDES FROM READS YOU WANT>
SAMPLE=${SMPL}_${NUCL}nt

printf "\nExtract UMIs plus first ${NUCL} nucleotides from read and make a UMI fastq\n\n"
printf "Start pipeline log\n"

start=$SECONDS
python3 $SCRIPTS/makeUmiFastq.py $FASTQ1 $OUTPUT_DIR/${SAMPLE}.umis.fastq $NUCL
end=$SECONDS
printf "\nDuration: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nMake unmapped BAM file ${SAMPLE}.unmapped.bam from the FASTQs\n\n"

start=$SECONDS
java -jar /hpc/cog_bioinf/common_scripts/picard-tools-1.98/FastqToSam.jar \
    F1=$FASTQ1 F2=$FASTQ2 O=$OUTPUT_DIR/${SAMPLE}.unmapped.bam SM=$SAMPLE TMP_DIR=$TMPDIR
end=$SECONDS
printf "\nDuration of Make unmapped BAM: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nAnnotate BAM file with UMIs from ${SAMPLE}.umis.fastq\n\n"

start=$SECONDS
java -jar $FGBIO --tmp-dir=$TMPDIR AnnotateBamWithUmis \
    --input=$OUTPUT_DIR/${SAMPLE}.unmapped.bam --output=$OUTPUT_DIR/${SAMPLE}.unmapped.UMIS.bam \
    --fastq=$OUTPUT_DIR/${SAMPLE}.umis.fastq
end=$SECONDS
printf "\nDuration of Annotate BAM with UMIs: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nGroup annotated BAM by UMIs\n\n"

start=$SECONDS
python3 $SCRIPTS/groupReadsByUmi.py $OUTPUT_DIR/${SAMPLE}.unmapped.UMIS.bam $OUTPUT_DIR/${SAMPLE}.unmapped.grouped.bam
end=$SECONDS
printf "\nDuration of Group Reads: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nSort grouped BAM on TemplateCoordinate\n\n"

start=$SECONDS
java -jar $FGBIO --tmp-dir=$TMPDIR SortBam \
    --input=${OUTPUT_DIR}/${SAMPLE}.unmapped.grouped.bam \
    --output=${OUTPUT_DIR}/${SAMPLE}.unmapped.grouped.sorted.bam \
    --sort-order=TemplateCoordinate
end=$SECONDS
printf "\nDuration of Sort Grouped BAM: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nCall molecular consensus reads of sorted and grouped BAM\n\n"

start=$SECONDS
java -Xmx4g -jar $FGBIO --tmp-dir=$TMPDIR \
     CallMolecularConsensusReads \
    --input=${OUTPUT_DIR}/${SAMPLE}.unmapped.grouped.sorted.bam \
    --output=${OUTPUT_DIR}/${SAMPLE}.consensus.unmapped.bam \
    --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=30 \
    --min-reads 1
end=$SECONDS
printf "\nDuration of Call Consensus: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nFilter molecular consensus reads with minimal depth of 3\n\n"

start=$SECONDS
java -Xmx4g -jar $FGBIO --tmp-dir=$TMPDIR \
    FilterConsensusReads \
    --input=${OUTPUT_DIR}/${SAMPLE}.consensus.unmapped.bam \
    --output=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.unmapped.bam \
    --ref=$REF_GENOME \
    --reverse-per-base-tags=true \
    -M 1 \
    -N 40 \
    -e 0.1 \
    -n 0.1
end=$SECONDS
printf "\nDuration of Filter Consensus: $((end-start)) seconds." >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nSort filtered consensus reads on read name\n"

start=$SECONDS
java -Xmx4g -jar $PICARD SortSam \
     I=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.unmapped.bam \
     O=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.unmapped.sorted.bam \
     SORT_ORDER=queryname TMP_DIR=$TMPDIR
end=$SECONDS
printf "\nDuration of Sort Filtered Consensus: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log

printf "\nGenerate a mapped BAM from the filtered, consensus reads\n\n"

start=$SECONDS
java -Xmx4g -jar $PICARD SamToFastq \
    I=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.unmapped.sorted.bam \
    F=/dev/stdout INTERLEAVE=true TMP_DIR=$TMPDIR \
    | bwa mem -p -M -t 8 $REF_GENOME /dev/stdin \
    | java -Xmx4g -jar ${PICARD} MergeBamAlignment \
        UNMAPPED=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.unmapped.sorted.bam \
        ALIGNED=/dev/stdin O=${OUTPUT_DIR}/${SAMPLE}.consensus.filtered.mapped.bam \
        R=${REF_GENOME} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
        ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$TMPDIR
end=$SECONDS
printf "\nDuration of Mapping: $((end-start)) seconds.\n" >> $OUTPUT_DIR/${SAMPLE}.time.log
printf "\nEnd of pipeline" >> $OUTPUT_DIR/${SAMPLE}.time.log

end1=$SECONDS

printf "\nDuration of whole pipeline is $((end-start)) seconds." >> $OUTPUT/${SAMPLE}.time.log
