# MBV Pipeline
# Jack Humphrey 
import glob
import pandas as pd
import os
import pathlib
R_VERSION = "R/4.0.3"
configfile: 'config.yaml'

dataCode = config["dataCode"]
VCF = config["VCF"]
bamFolder = config["bamFolder"]
bamSuffix = config["bamSuffix"]
out_folder = config["out_folder"]
BAM_SAMPLES = [os.path.basename(x).replace(bamSuffix, "") for x in glob.glob(bamFolder + "/*" + bamSuffix)]
dataCode = dataCode + "_mbv"

#bamstats_output = "bamstats/{sample}.bamstat.txt"
final_output = out_folder + dataCode + "_summary.txt"

print(" * MBV pipeline *")
print(" Jack Humphrey 2019-2023")
print(" * Data code is : %s " % dataCode)

#pathlib.Path(out_folder + "bamstats").mkdir(parents=True, exist_ok=True)

print(BAM_SAMPLES)
print(len(BAM_SAMPLES) )
rule all:
    input:
        #expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.{group_by}.combined_covariates.txt", PEER_N = PEER_values, group_by = group_by_values)
        final_output

rule indexBam:
    input:
        bamFolder + "{sample}.bam"
    output:
        bamFolder + "{sample}.bam.bai"
    shell:
        "ml samtools/1.9;"
        "samtools index {input}"

# run match BAM to variants
rule matchBAM2VCF:
    input:
        bam = bamFolder + "{sample}.bam",
        bai = bamFolder + "{sample}.bam.bai",
        vcf = VCF
    output:
        out_folder + "{sample}.bamstat.txt"
    shell:
        "ml qtltools/1.2;"
        "QTLtools mbv --bam {input.bam} --vcf {input.vcf} --filter-mapping-quality 150 --out {output}"

#rule summariseResults:
#    input:
#        files = expand(out_folder + "{sample}.bamstat.txt", sample = BAM_SAMPLES)
#    output:
#        final_output
#    shell:
#        "set +o pipefail;"
#        "for i in {input.files};"
#        "do cat $i | sort -k9nr,10nr | head -1 | awk -v i=$i \'{{print i, $0}}\'  ;"
#        "done > {output};"

rule collateMBV:
    input:
         files = expand(out_folder + "{sample}.bamstat.txt", sample = BAM_SAMPLES)
    output:
        final_output
    shell:
        "ml {R_VERSION};"
        "Rscript scripts/collate_mbv.R --outFolder {out_folder} --outFile {final_output}"
