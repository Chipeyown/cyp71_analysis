bowtie -n 0 -v 0 -k 1 /BioII/qiyijun/ltf_20181126/DBs/bowtie_index/ath/ath10  -p 12 -q fastq_processed/$pre.trim.fq -S >$pre.sam 2>$pre.mapped.log
