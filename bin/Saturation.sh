#!/bin/bash
for i in 1 2 3 4 5 6 7 8 9
do
        samtools view -@40 -s 100.${i} ./${2} -b -o ${1}_${i}0%.bam
        featureCounts -T 40 -f -t gene -p -B -C -g gene_id -a ${3} -o ${1}_counts.txt ${1}_${i}0%.bam
        count=($(awk '{print $7}' ${1}_counts.txt))
        num=${#count[@]}                   
        nonzero=0
        for((j=0;j<$num;j++));do
                if [ ${count[$j]} != 0 ];then
                        let nonzero++
                fi
        done
        echo "${1}"     "${1}_${i}%"    $nonzero >> ${1}.Gene_Number.txt
        echo "${1}_${i}%"       `samtools view -@40 -c ${1}_${i}0%.bam` >> ${1}.Reads_Number.txt


done
