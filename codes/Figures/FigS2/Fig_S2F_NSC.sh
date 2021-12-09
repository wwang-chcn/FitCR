#! /bin/bash

cd ~/fitCUTRUN/TF/encode_data_quality/bin

## =======================
### ATF1 by sample (merge replicate)
### =======================

cd ~/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample
mkdir ATF1_FitCR_1e5 && cd ATF1_FitCR_1e5
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ATF1_C_1e5_rep1.bam
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ATF1_C_1e5_rep2.bam
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ATF1_C_1e5_rep3.bam
cd ../
cd ATF1_FitCR_1e5
# samtools merge -f ATF1_1e5.bam ATF1_C_1e5_rep1.bam ATF1_C_1e5_rep2.bam ATF1_C_1e5_rep3.bam &
cd ../

mkdir ATF1_FlagCR_1e5 && cd ATF1_FlagCR_1e5
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ATF1_1e5_rep1.bam ATF1_FlagCR_1e5_rep1.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ATF1_1e5_rep2.bam ATF1_FlagCR_1e5_rep2.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ATF1_1e5_rep3.bam ATF1_FlagCR_1e5_rep3.bam
cd ../
cd ATF1_FlagCR_1e5
samtools merge -f ATF1_FlagCR_1e5.bam ATF1_FlagCR_1e5_rep1.bam ATF1_FlagCR_1e5_rep2.bam ATF1_FlagCR_1e5_rep3.bam &
cd ../

mkdir ATF1_wtCR_1e5 && cd ATF1_wtCR_1e5
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ATF1_1e5_rep1.bam ATF1_wtCR_1e5_rep1.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ATF1_1e5_rep2.bam ATF1_wtCR_1e5_rep2.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ATF1_1e5_rep3.bam ATF1_wtCR_1e5_rep3.bam
cd ../
cd ATF1_wtCR_1e5
samtools merge -f ATF1_wtCR_1e5.bam ATF1_wtCR_1e5_rep1.bam ATF1_wtCR_1e5_rep2.bam ATF1_wtCR_1e5_rep3.bam &
cd ../

mkdir ATF1_ChIP && cd ATF1_ChIP
ln -s ~/fitCUTRUN/TF/encode_data_quality/0_raw_data/atf1_K562_rep1.fastq.gz .
ln -s ~/fitCUTRUN/TF/encode_data_quality/0_raw_data/atf1_K562_rep2.fastq.gz .
cd ../
cd ATF1_ChIP
cat atf1_K562_rep1.fastq.gz atf1_K562_rep2.fastq.gz > ATF1_ChIP.fastq.gz
cd ../


## =======================
### ELF1 by sample (merge replicate)
### =======================

cd ~/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample
mkdir ELF1_FitCR_1e5 && cd ELF1_FitCR_1e5
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ELF1_C_1e5_rep1.bam
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ELF1_C_1e5_rep2.bam
ln -s /mnt/Storage/home/wangwen/project/zebrafish_FitCR/202103/1_mapping/ELF1_C_1e5_rep3.bam
cd ../
cd ELF1_FitCR_1e5
samtools merge -f ELF1_1e5.bam ELF1_C_1e5_rep1.bam ELF1_C_1e5_rep2.bam ELF1_C_1e5_rep3.bam &
cd ../

mkdir ELF1_FlagCR_1e5 && cd ELF1_FlagCR_1e5
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ELF1_1e5_rep1.bam ELF1_FlagCR_1e5_rep1.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ELF1_1e5_rep2.bam ELF1_FlagCR_1e5_rep2.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-mFlag_ELF1_1e5_rep3.bam ELF1_FlagCR_1e5_rep3.bam
cd ../
cd ELF1_FlagCR_1e5
samtools merge -f ELF1_FlagCR_1e5.bam ELF1_FlagCR_1e5_rep1.bam ELF1_FlagCR_1e5_rep2.bam ELF1_FlagCR_1e5_rep3.bam &
cd ../

mkdir ELF1_wtCR_1e5 && cd ELF1_wtCR_1e5
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ELF1_1e5_rep1.bam ELF1_wtCR_1e5_rep1.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ELF1_1e5_rep2.bam ELF1_wtCR_1e5_rep2.bam
ln -s /mnt/Storage/home/wangyiman/fitCUTRUN/TF/2021*/1_mapping/CR-wt_ELF1_1e5_rep3.bam ELF1_wtCR_1e5_rep3.bam
cd ../
cd ELF1_wtCR_1e5
samtools merge -f ELF1_wtCR_1e5.bam ELF1_wtCR_1e5_rep1.bam ELF1_wtCR_1e5_rep2.bam ELF1_wtCR_1e5_rep3.bam &
cd ../

mkdir ELF1_ChIP && cd ELF1_ChIP
ln -s ~/fitCUTRUN/TF/encode_data_quality/0_raw_data/elf1_K562_rep3.fastq.gz .
ln -s ~/fitCUTRUN/TF/encode_data_quality/0_raw_data/elf1_K562_rep4.fastq.gz .
cd ../
cd ELF1_ChIP
cat elf1_K562_rep3.fastq.gz elf1_K562_rep4.fastq.gz > ELF1_ChIP.fastq.gz &
cd ../


cd ~/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample


bash strand_cross_corr.sh ATF1_1e5.bam ATF1_FitCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ATF1_FitCR_1e5 false
bash strand_cross_corr.sh ELF1_1e5.bam ELF1_FitCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ELF1_FitCR_1e5 false

bash strand_cross_corr.sh ATF1_FlagCR_1e5.bam ATF1_FlagCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ATF1_FlagCR_1e5 false
bash strand_cross_corr.sh ELF1_FlagCR_1e5.bam ELF1_FlagCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ELF1_FlagCR_1e5 false

bash strand_cross_corr.sh ATF1_wtCR_1e5.bam ATF1_wtCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ATF1_wtCR_1e5 false
bash strand_cross_corr.sh ELF1_wtCR_1e5.bam ELF1_wtCR /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ELF1_wtCR_1e5 false

bash strand_cross_corr.sh ATF1_ChIP.fastq.gz ATF1_ChIP /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ATF1_ChIP false
bash strand_cross_corr.sh ELF1_ChIP.fastq.gz ELF1_ChIP /mnt/Storage/home/wangyiman/source/bySpecies/hg38/bowtie2/hg38_main 5 /mnt/Storage/home/wangyiman/fitCUTRUN/revision_1/analysis/TF/comparison/NSC_bySample/ELF1_ChIP false


