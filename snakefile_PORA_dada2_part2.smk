## env env_variants
## env dada2

##snakemake --snakefile scripts/PORA_dada2.smk --configfile configs/config.yml --cores 20

import pandas as pd
import glob
import sys
from pathlib import Path
from Bio import SeqIO


## --------------------------------------------------------------------------------
## global parameters from config file
PREFIX=config["folder"]
RUN=config["runid"]

## --------------------------------------------------------------------------------
## targets - rules to call
remove_primers_complete=expand(PREFIX+"/Resources/data/primers_removed/" + RUN + "/{sample}_rmprimer_R1.fastq.gz", sample=config["samplesR1"])
process=expand(PREFIX+ "{files}", files=["/Results/data/seqtabs/" + RUN + "/seqtab.rds",
        "/Results/plots/" + RUN + "/Read1_quality.pdf",
        "/Results/plots/" + RUN + "/Read2_quality.pdf",
        "/Results/plots/" + RUN + "/filter_read1_quality.pdf",
        "/Results/plots/" + RUN + "/filter_read2_quality.pdf",
        "/Results/plots/" + RUN + "/Read1_errorRate.pdf",
        "/Results/plots/" + RUN + "/Read2_errorRate.pdf",
        "/Results/tables/" + RUN + "/reads_overview.tsv"])
#process=expand(PREFIX+ "/plots/" + RUN + "/{sample}_read1_quality.pdf", sample=config["samplesR1"])
#process=expand(PREFIX+ "/plots/" + RUN + "/Read1_quality.pdf")



rule all:
    input:
        remove_primers_complete + process

## --------------------------------------------------------------------------------
## rules

wildcard_constraints:
    sample=r"[\w\d-]+"

rule dada2:
    input:
        r1=expand(PREFIX+"/Resources/data/primers_removed/" + RUN + "/{sample}_rmprimer_R1.fastq.gz", sample=config["samplesR1"]), 
        r2=expand(PREFIX+"/Resources/data/primers_removed/" + RUN + "/{sample}_rmprimer_R2.fastq.gz", sample=config["samplesR1"]), 
    output:
        PREFIX+"/Results/data/seqtabs/" + RUN + "/seqtab.rds", 
        PREFIX+"/Results/plots/" + RUN + "/Read1_quality.pdf",
        PREFIX+"/Results/plots/" + RUN + "/Read2_quality.pdf",
        PREFIX+"/Results/plots/" + RUN + "/filter_read1_quality.pdf",
        PREFIX+"/Results/plots/" + RUN + "/filter_read2_quality.pdf",
        PREFIX+"/Results/plots/" + RUN + "/Read1_errorRate.pdf",
        PREFIX+"/Results/plots/" + RUN + "/Read2_errorRate.pdf",
        PREFIX+"/Results/tables/" + RUN + "/reads_overview.tsv",
    threads: 
        1
    conda:
        'FBI_amplicons'
    params:
        prefix=PREFIX,
        run=RUN,
        trunc1=config["trunc1"],
        trunc2=config["trunc2"],
    shell:
        """ 
        dada2_auto.R --r1 $(realpath --relative-to="$PWD" {input.r1} | tr '\\n' ',') --r2 $(realpath --relative-to="$PWD" {input.r2} | tr '\\n' ',') --prefix {params.prefix}/Results --run {params.run} --trunc1 {params.trunc1} --trunc2 {params.trunc2}
        """
