# MBV_Pipeline

MBV is part of the QTLtools pakcage. It ensures good match between the sequence and genotype data. This is useful to detect sample mislabeling, contamination or PCR amplification biases.

This pipeline parallelizes MBV on the cluster using Snakemake.

To Run:

1. First edit the `config.yaml` file asnecessary. The configurations are as follows:
    - `bamFolder`: A folder containing all the aligned RNA-seq data to be matched. These should be `bam` files.
    - `bamSuffix`: The suffix of your input bam files. This will probably be ".bam" so just leave this one alone.
    - `out_folder`: The output directory where you want your final output to go. "output/" works or you can come up with a more exciting name.
    - `dataCode`: The name of your data to be used for identification purposes later on. Pick something non-generic, e.g. ""
    - `VCF`: The path to the genotype data to matched with the above RNA-seq data. This should be one `vcf` file containing all chromosomes.

2. Run Snakemake
first activate the snakemake environment
```
conda activate snakemake
```

then perform a dry run to make sure all the files are in place and to check that the output will be correct (Snakefile assumes config.yaml is the correct config file, however, a new one can be specified with the option `--configfile <pick_a_file.whatever>`)
```
snakemake -s Snakefile -npr
```

If everything looks good, you can now run the pipeline (this assumes you have your environment set up properly - see https://github.com/RajLabMSSM/RNA-pipelines/tree/master/snakejob on how to do that)
```
snakejob_HPC -s Snakemake
```


And that should really be it!

The final output should consist of a summary table of all the match information called < `dataCode` >_mbv_summary.txt
