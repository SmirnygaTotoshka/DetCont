#!/bin/bash

HOME_DIR="/media/HEAP-EPI/stotoshka/DetCont"
SUB="epi"

cd $HOME_DIR

for dir in $(find $HOME_DIR -maxdepth 1 -type d)
do
    if grep -q "$SUB" <<< "$dir"; then
        cd $dir
        bam_file=$(ls | grep ".bam$")
        fastq_file=$(basename -s .bam $bam_file)".fastq"
        echo $dir
        bedtools bamtofastq -i $bam_file -fq $fastq_file
    fi
done
