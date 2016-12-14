from validation import perform_sims, run_snppipe, perform_analyses


dirstub = "validation/lm_full/run"
config = "validation/lm/Lm.config"
refloc = "validation/lm/CFSAN023463.fasta"
num = 5

#perform_sims(config, dirstub, num)
#run_snppipe(dirstub, refloc, num)
perform_analyses(dirstub, num)


