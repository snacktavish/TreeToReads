
TTR_OUT=/home/ejmctavish/Documents/FDA/TreetoReads/test_out
DATA=cfsan_snppipe_input
REF=/home/ejmctavish/Documents/FDA/TreetoReads/example/barref.fasta

cd $TTR_OUT

mkdir $DATA

for fi in `ls ${TTR_OUT}/*_1.fq | cut -d"," -f2`
    do 
        fullp=${fi%_*}
        filename=${fullp##*/}
        echo $filename
        mkdir $DATA/$filename
        cp ${fullp}_1.fq $DATA/$filename
        cp ${fullp}_2.fq $DATA/$filename
        
    done

    
run_snp_pipeline.sh -s $DATA $REF