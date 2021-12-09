#!/bin/bash

cd $word_dir

function sampling {
    n=`cat ${rep*_reads_bed} | wc -l `
    ratio=`bc -l <<< "10000000 / ${n}"`
    cat ${rep*_reads_bed} | awk -v ratio=${ratio} 'BEGIN{srand(1006)} {if(rand()<ratio) print $0}' > ${sampling_down_reads_bed}
}

function peak_calling {
    macs2 callpeak -f BEDPE -t ${sampling_down_reads_bed} -n ${TF_stage_name} -g hs --outdir $work_dir -q 0.01 --keep-dup all 2>&1 >>/dev/null | tee ${MACS_log}
}

function piling_up {
    fragment_length=`awk 'BEGIN{s=0} {s+=$3-$2} END{print s/NR}' ${sampling_down_reads_bed}`
    ShiftPairEnd.sh ${sampling_down_reads_bed} ${fragment_length}
    n=`wc -l ${sampling_down_reads_bed_shift} | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${sampling_down_reads_bed_shift} -g danRer11_2_main.chrom.sizes > ${sampling_down_reads_bed_shift_bdg} && \
    bdg2bw.sh ${sampling_down_reads_bed_shift_bdg} danRer11_2_main.chrom.sizes ${TF_stage_name} && \
    rm ${sampling_down_reads_bed_shift} ${sampling_down_reads_bed_shift_bdg}
}

qvalue_cutoff=10
fold_cutoff=5

function filter {
    awk -v fold_cutoff=${fold_cutoff} -v qvalue_cutoff=${qvalue_cutoff} 'NR>21{if($8>fold_cutoff && $9>qvalue_cutoff) printf "%s\t%d\t%d\t%s\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n", $1, $2-1, $3, $10, $9*10, $8, $7, $9, $5-$2}' ${TF_stage_name_peaks_xls} > ${TF_stage_name_filtered_peaks}
    awk -v fold_cutoff=${fold_cutoff} -v qvalue_cutoff=${qvalue_cutoff} 'NR>21{if($8>fold_cutoff && $9>qvalue_cutoff) printf "%s\t%d\t%d\t%s\t%.5f\n", $1, $5-1, $5, $10, $9}' ${TF_stage_name_peaks_xls} > ${TF_stage_name_filtered_summits}
}

function exclude_IgG {
    bedtools intersect -v -a ${TF_stage_name_filtered_peaks} -b ${IgG_filtered_peaks} > ${TF_stage_name_filtered_peaks_exclude_IgG} 

}

sampling && peak_calling && piling_up && filtering

if $tech == 'FitCUT&RUN' | $tech == 'CUT&RUN' | $tech == 'FLAG-CUT&RUN';then
    exclude_IgG
fi
