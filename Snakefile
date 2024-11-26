# MBV Pipeline
# Jack Humphrey 
import glob
import pandas as pd
import os
import pathlib
R_VERSION = "R/4.0.3"
#somalier="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/somalier"
somalier="/sc/arion/projects/ad-omics/data/software/somalier_0.2.19/somalier"
configfile: 'config.yaml'

dataCode = config["dataCode"]
VCF = config["VCF"]
bamFolder = config["bamFolder"]
bamSuffix = config["bamSuffix"]
out_folder = config["out_folder"]
BAM_SAMPLES = [os.path.basename(x).replace(bamSuffix, "") for x in glob.glob(bamFolder + "/*" + bamSuffix)]
dataCode = dataCode + "_mbv"

som_folder = "somalier/"

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
        final_output #,
        #som_folder + dataCode + ".somalier-ancestry.html"

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
        "ml qtltools/1.3;"
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

## SOMALIER RULES

rule somalier_bam:
    input:
        bam = bamFolder + "{sample}.bam",
        bai = bamFolder + "{sample}.bam.bai"
    params:
        sites="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/sites/sites.hg38.vcf.gz",
        ref="/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa",
        ancestry_labels="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/ancestry-labels-1kg.tsv",
        labeled_samples="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/1kg-somalier/"
    output:
        som_folder + "{sample}.somalier"
    shell:
        "export SOMALIER_SAMPLE_NAME={wildcards.sample};"
        "{somalier} extract -d {som_folder} --sites {params.sites} -f {params.ref} {input.bam}"

# not using for now
rule somalier_vcf:
    input:
        vcf = VCF
    params:
        sites="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/sites/sites.hg38.vcf.gz",
        ref="/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa",
        ancestry_labels="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/ancestry-labels-1kg.tsv",
        labeled_samples="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/1kg-somalier/"
    output:
        expand(som_folder + "{sample}.somalier", sample = BAM_SAMPLES),
        som_folder + "somalier_complete.txt"
    shell:
        "{somalier} extract -d {som_folder} --sites {params.sites} -f {params.ref} {input.vcf};"
        "touch {som_folder}/somalier_complete.txt"



rule relatedness_ancestry:
    input:
        expand(som_folder + "{sample}.somalier", sample = BAM_SAMPLES)
        #som_folder +  "somalier_complete.txt"
    output:
        som_folder + dataCode + ".somalier-ancestry.html",
        som_folder + dataCode + ".somalier.relatedness.html",
        som_folder + dataCode + ".somalier-ancestry.tsv"
    params:
        sites="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/sites/sites.hg38.vcf.gz",
        ref="/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa",
        ancestry_labels="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/ancestry-labels-1kg.tsv",
        labeled_samples="/sc/arion/projects/ad-omics/data/software/somalier_0.2.12/ancestry/1kg-somalier/"
    shell:
        "{somalier} relate {som_folder}*.somalier;"
        "{somalier} ancestry --labels {params.ancestry_labels} {params.labeled_samples}*.somalier ++ {som_folder}*.somalier;"
        "rename somalier {dataCode}.somalier somalier* ;"
        "mv {dataCode}.somalier* {som_folder};"
         #      mv {dataCode}.somalier.html {dataCode}.somalier.relatedness.html; \
          #     mv {dataCode}.somalier-ancestry.somalier-ancestry.tsv {dataCode}.somalier.ancestry.tsv; \
           #    mv {dataCode}.somalier-ancestry.somalier-ancestry.html {dataCode}.somalier.ancestry.html; \
            # ')
