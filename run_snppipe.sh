#!/bin/sh
#Make sure FDA/snp-pipeline/bin is in the path

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
        cp ${fullp}_1.fq $TTR_OUT/$DATA/$filename
        cp ${fullp}_2.fq $TTR_OUT/$DATA/$filename
        
    done

    
run_snp_pipeline.sh -s $TTR_OUT/$DATA -o $TTR_OUT/$DATA $REF

cp example/garli.conf $TTR_OUT/$DATA

cd $TTR_OUT/$DATA

Garli 