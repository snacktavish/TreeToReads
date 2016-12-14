from validation import perform_sims, run_snppipe, perform_analyses
import os

dirstub = "validation/barref_full/run"
config = "validation/tree1.config"
refloc =  "example/barref.fasta"
altrefloc = 'validation/LT2.fasta'
num = 5

#perform_sims(config, dirstub, num)

#for i in range(num):
#        i+=1
#        diri = "{}{}".format(dirstub,i)
#        os.system("mkdir {}/altref".format(diri))
#        os.system("cp -r {}/fastq {}/altref/fastq".format(diri, diri))


#run_snppipe(dirstub, refloc, num)
inf_dict = perform_analyses(dirstub, num)



#alt refreence genome
'''
for i in range(num):
    i+=1
    diri = "{}{}".format(dirstub,i)
    os.system("run_snp_pipeline.sh -f -s {}/altref/fastq -o {}/altref {}".format(diri, diri, altrefloc)) 
'''

