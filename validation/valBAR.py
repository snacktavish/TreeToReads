from validation import perform_sims, perform_analyses


dirstub = "validation/barref_full/run"
config = "validation/tree1.config"
#perform_sims(config, dirstub, "example/barref.fasta", 10)
perform_analyses(dirstub, 10)
