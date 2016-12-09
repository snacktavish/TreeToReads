from validation import perform_sims, perform_analyses


dirstub = "validation/lm_full/run"
config = "validation/lm/Lm.config"
perform_sims(config, dirstub, "validation/lm/CFSAN023463.fasta", 10)
perform_analyses(dirstub, 10)


