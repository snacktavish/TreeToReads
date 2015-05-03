#!/bin/sh
#Make sure FDA/snp-pipeline/bin is in the path

TTR_OUT=$1
REF=$2
DATA=cfsan_snppipe_input_test

#CHECK ARGS!
cd $TTR_OUT
mkdir -p $DATA #How to make this a try statement?

for fi in `ls *_1.fq | cut -d"," -f2`
    do 
        fullp=${fi%_*}
        filename=${fullp##*/}
        echo $filename
        mkdir -p $DATA/$filename
        cp ${fullp}_1.fq $DATA/$filename
        cp ${fullp}_2.fq $DATA/$filename
        
    done

    
run_snp_pipeline.sh -s $DATA -o nee $REF

cp ../example/garli.conf $DATA

cd $DATA

Garli 