#!/bin/sh


TTR_OUT=$1
REF=$2
DATA=cfsan_snppipe_input_test

#CHECK ARGS!

#mkdir $DATA #How to make this a try statement?

for fi in `ls ${TTR_OUT}/*_1.fq | cut -d"," -f2`
    do 
        fullp=${fi%_*}
        filename=${fullp##*/}
        echo $filename
        mkdir $DATA/$filename
        cp ${fullp}_1.fq $DATA/$filename
        cp ${fullp}_2.fq $DATA/$filename
        
    done

    
#run_snp_pipeline.sh -s $DATA $REF # need to remeber how to install SNP pipeline...