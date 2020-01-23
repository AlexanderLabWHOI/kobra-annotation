configfile: "config.yaml"

import pandas as pd
import numpy as np
import os
import io

KEGG_DB = config['kegg_db']
KEGG_DIR = config['kegg_dir']
BLASTPREFIX = os.path.join(KEGG_DIR, 'blastdb', os.path.split(KEGG_DB)[-1]).strip('.pep')
BLASTDB = expand(BLASTPREFIX+'.00.{ext}', ext =["phr", "pin", "psq"])

ID, = glob_wildcards("genomes/{id}.fa")

print(ID)
#---RULES---#

rule all:
    input: BLASTDB, expand('output/blast/{id}.blast.tab', id=ID)

rule build_blast_db:
    input: KEGG_DB
    output: BLASTDB
    params: 
        prefix= BLASTPREFIX
    conda: 
        'envs/blast.yaml'
    shell:
        """
        makeblastdb -dbtype prot -in {input} -out {params.prefix}
        """

rule run_blastx:
    input: 'genomes/{id}.fa'
    output: 'output/blast/{id}.blast.tab'
    params:
        db = BLASTPREFIX
    conda:
        'envs/blast.yaml'

    shell:
        """
        blastx -db {params.db} -query {input} -outfmt 6 -evalue 1e-5
        """
