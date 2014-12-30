
TTR_OUT=/home/ejmctavish/Documents/FDA/TreetoReads/test_out
PIPEFI=cfsan_snppipe_input
REF=/home/ejmctavish/Documents/FDA/TreetoReads/example/barref.fasta

cd $TTR_OUT

mkdir $PIPEFI

for fi in `ls ${TTR_OUT}/*_1.fq | cut -d"," -f2`
    do 
        fullp=${fi%_*}
        filename=${fullp##*/}
        echo $filename
        mkdir $PIPEFI/$filename
        cp ${fullp}_1.fq $PIPEFI/$filename
        cp ${fullp}_2.fq $PIPEFI/$filename
        
    done

    
run_snp_pipeline.sh -s $PIPEFI $REF