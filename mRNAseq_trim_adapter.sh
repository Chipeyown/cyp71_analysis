#ungzip
gunzip $file1
gunzip $file2
touch $pre.inputlist
echo $file1>>$pre.inputlist
echo $file2>>$pre.inputlist
fastuniq -i $pre.inputlist -t q -o $file1 -p $file2 -c 0
trim_galore --clip_R1 9 --clip_R2 9 --stringency 2 -o ../fastq_processed --paired $file1 $file2