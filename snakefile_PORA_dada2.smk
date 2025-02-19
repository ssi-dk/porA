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

rule all:
    input:
        remove_primers_complete 

## --------------------------------------------------------------------------------
## rules

wildcard_constraints:
    sample=r"[\w\d-]+"


rule remove_primers:
    input:
        r1=lambda wildcards: config["samplesR1"][wildcards.sample],
        r2=lambda wildcards: config["samplesR2"][wildcards.sample],
    output:
        r1=PREFIX+"/Resources/data/primers_removed/" + RUN + "/{sample}_rmprimer_R1.fastq.gz",
        r2=PREFIX+"/Resources/data/primers_removed/" + RUN + "/{sample}_rmprimer_R2.fastq.gz"
    log:
        PREFIX+"/Resources/data/cutadapt_logs/" + RUN + "/{sample}_rmprimer.log"
    threads: 
        4
    params:
        part=config["partition"],
        fwd="CCACAATTATGGTTAGCTTA",
        rev="TGAGAAGTTAAGTTTTGGAGAG"

    shell:
        """
        srun -p {params.part} --mem 10000 -c {threads} cutadapt -g {params.fwd} -G {params.rev} -e 0.08 -j {threads} -m 30 --discard-untrimmed --pair-adapters --pair-filter both -o {output.r1} -p {output.r2} {input.r1} {input.r2}&> {log} 
        """


