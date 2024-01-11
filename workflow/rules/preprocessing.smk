rules preprocess_PARIS:
    shell:
        remove 5' and 3' adapters using trimmomatic-0.32.jar SE -threads 16 -phred33
removed PCR duplicates using splitFastqLibrary (in the icSHAPE pipeline)