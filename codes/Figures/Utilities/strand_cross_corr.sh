#!/bin/bash

### bash strand_cross_corr.sh <FASTQ_R1|BAM_PE> <PREFIX> <BWT2_IDX> <NTH_BWT2> <WORKDIR> <OCR_FLAG>

PREFIX=$2
BWT2_IDX=$3
NTH_BWT2=$4 ### number of threads of bowtie2
WORKDIR=$5
OCR_FLAG=$6

cd $WORKDIR

case ${1##*.} in
    "gz")
        echo "SE FASTQ mode :"
        FASTQ_R1=$1
        # =================================
        ### Trim R1 fastq to 50bp
        # =================================
        TRIMMED_FASTQ_R1=trimmed_${FASTQ_R1}
        RAW_READ_LEN=$(zcat $FASTQ_R1 | head -n 100 | awk 'BEGIN{b=0}NR%4==2{a=length;b=b+a}END{print b/NR*4}')
        a=${RAW_READ_LEN%%.*}
        b=50
        TRIMMED_RAW_READ_LEN=$(( a < b ? a : b ))

        [ ! -f $TRIMMED_FASTQ_R1 ] && python /mnt/Storage/home/wangyiman/software/chip-seq-pipeline2/src/trimfastq.py $FASTQ_R1 $TRIMMED_RAW_READ_LEN | gzip -nc >  $TRIMMED_FASTQ_R1


        # =================================
        ### Align $TRIMMED_FASTQ_R1 (not paired) with bowtie2 (step 1a SE) and use it for filtering step (1b) and then get $FILT_BAM_FILE (not the deduped $FINAL_BAM_FILE), which is filtered but not deduped.
        # =================================
        RAW_BAM_FILE=${PREFIX}.bam
        log=${PREFIX}.align.log
        flagstat_qc=${PREFIX}.flagstat.qc
        bowtie2 --mm -x $BWT2_IDX --threads $NTH_BWT2 -U $TRIMMED_FASTQ_R1 2> $log | /mnt/Storage/home/wangyiman/anaconda3/envs/encode-chip-seq-pipeline/bin/samtools view -Su /dev/stdin | /mnt/Storage/home/wangyiman/anaconda3/envs/encode-chip-seq-pipeline/bin/samtools sort - -o $RAW_BAM_FILE
        /mnt/Storage/home/wangyiman/anaconda3/envs/encode-chip-seq-pipeline/bin/samtools sort -n --threads 10 ${RAW_BAM_FILE} -O SAM | SAMstats --sorted_sam_file - --outf ${flagstat_qc}

        FILT_BAM_PREFIX="${PREFIX}.filt.srt"
        FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
        MAPQ_THRESH=30
        /mnt/Storage/home/wangyiman/anaconda3/envs/encode-chip-seq-pipeline/bin/samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} -o ${FILT_BAM_FILE}

        echo "Align is done!"

        TA_FILE=${PREFIX}.SE.tagAlign.gz
        bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > ${TA_FILE}
        ;;

    "bam")
        echo "PE BAM mode :"
        BAM_PE=$1
        TA_FILE=${PREFIX}.SE.tagAlign.gz
        case $OCR_FLAG in
            "true")
                bamToBed -bedpe -mate1 -i $BAM_PE | awk '$1 !~ /_/{f_len=0; if($2<$5) {f_len=$6-$2} else {f_len=$3-$5}; if (f_len<=120) {if ($9 == "+"){if ($3-$2<=50){print $1"\t"$2"\t"$3"\t""N""\t""1000""\t"$9} else {print $1"\t"$2"\t"$2+50"\t""N""\t""1000""\t"$9}} else {if ($3-$2<=50){print $1"\t"$2"\t"$3"\t""N""\t""1000""\t"$9} else {print $1"\t"$3-50"\t"$3"\t""N""\t""1000""\t"$9}}}}' > ${TA_FILE%.gz}
                ;;
            "false")
                bamToBed -bedpe -mate1 -i $BAM_PE | awk '$1 !~ /_/{f_len=0; if($2<$5) {f_len=$6-$2} else {f_len=$3-$5}; if ($9 == "+"){if ($3-$2<=50){print $1"\t"$2"\t"$3"\t""N""\t""1000""\t"$9} else {print $1"\t"$2"\t"$2+50"\t""N""\t""1000""\t"$9}} else {if ($3-$2<=50){print $1"\t"$2"\t"$3"\t""N""\t""1000""\t"$9} else {print $1"\t"$3-50"\t"$3"\t""N""\t""1000""\t"$9}}}' > ${TA_FILE%.gz}
                ;;
        esac

        gzip -nc ${TA_FILE%.gz} > ${TA_FILE}
        ;;

esac



# =================================
# make tagAlign for filtered (but not deduped) BAM
# and subsample it for cross-correlation analysis
# ================================

NREADS=15000000
SUBSAMPLED_TA_FILE="${PREFIX}.filt.sample.$((NREADS / 1000000)).SE.tagAlign.gz"

zcat ${TA_FILE} | grep -v “chrM” | shuf -n ${NREADS} --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | gzip -nc > ${SUBSAMPLED_TA_FILE}

echo "Tag is done!"



# =================================
# Estimate read length from first 100 reads.
# =================================
zcat ${TA_FILE} > ${TA_FILE}.tmp
READ_LEN=$(head -n 100 ${TA_FILE}.tmp | awk 'function abs(v) {{return v < 0 ? -v : v}} BEGIN{{sum=0}} {{sum+=abs($3-$2)}} END{{print int(sum/NR)}}')


# =================================
# Determine exclusion range for fragment length estimation.
# Use a fixed lowerbound at -500.
# Upperbound EXCLUSION_RANGE_MAX is
#   TF ChIP-seq:  max(read_len + 10, 50)
#   Histone ChIP-seq:  max(read_len + 10, 100)
# lowerbound is fixed at 500 for both
# =================================
EXCLUSION_RANGE_MIN=-500
a=`bc -l <<< "${READ_LEN%%.*} + 10"`
echo Average read length is ${a} - 10.
b=50
EXCLUSION_RANGE_MAX=$(( a > b ? a : b ))
rm -f ${TA_FILE}.tmp


# =================================
### cross-correlation analysis
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
# =================================
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
CC_RDATA_FILE="${SUBSAMPLED_TA_FILE}.cc.RData"
NTHREADS=$NTH_BWT2

Rscript /mnt/Storage/home/wangyiman/software/phantompeakqualtools/run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -savd=${CC_RDATA_FILE} -out=${CC_SCORES_FILE} -x=${EXCLUSION_RANGE_MIN}:${EXCLUSION_RANGE_MAX} -s=-50:5:400 -rf 2>&1 >>/dev/null | tee ${SUBSAMPLED_TA_FILE}_cc.log

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}

echo "c-c is done!"



