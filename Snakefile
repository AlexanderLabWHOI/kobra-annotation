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
    input: BLASTDB, 
       # expand('output/blast/{id}.blastx.tab', id=ID),
        expand('output/blast/{id}.blastp.tab', id=ID)
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
    output: 'output/blast/{id}.blastx.tab'
    params:
        db = BLASTPREFIX
    conda:
        'envs/blast.yaml'

    shell:
        """
        blastx -db {params.db} -query {input} -out {output} -outfmt 6 -evalue 1e-5
        """

rule run_blastp: 
    input: 'genomes/{id}.faa'
    output: 'output/blast/{id}.blastp.tab'
    params:
        db = BLASTPREFIX
    conda:
        'envs/blast.yaml'
    shell:
        """
        blastp -db {params.db} -query {input} -out {output} -outfmt 6 -evalue 1e-5
        """


#rule kegg_annot:
#    input: 'output/blast/{id}.blastp.tab'
#    output: 'test'
#    conda:
#        'envs/kegg-annot.yaml'
#    shell:
#        """
#        keggannot_genes2ko
#        """
