configfile: "config.yaml"

import pandas as pd
import numpy as np
import os
import io

KEGG_DB = config['kegg_db']
KEGG_DIR = config['kegg_dir']
BLASTPREFIX = os.path.join(KEGG_DIR, 'blastdb', os.path.split(KEGG_DB)[-1]).strip('.pep')
BLASTDB = expand(BLASTPREFIX+'.00.{ext}', ext =["phr", "pin", "psq"])
DIAMOND_OUT = os.path.join(KEGG_DIR, 'diamond_db', os.path.split(KEGG_DB)[-1]).strip('.pep')
ID, = glob_wildcards("genomes/{id}.fa")

print(ID)

#---RULES---#

rule all:
    input: BLASTDB, 
       # expand('output/blast/{id}.blastx.tab', id=ID),
        expand('output/blast/{id}.blastp.tab', id=ID), expand('output/kegg_table/{id}.kegg_table.tsv', id=ID)

#rule build_blast_db:
#    input: KEGG_DB
#    output: BLASTDB
#    params: 
#        prefix= BLASTPREFIX
#    conda: 
#        'envs/blast.yaml'
#    shell:
#        """
#        makeblastdb -dbtype prot -in {input} -out {params.prefix}
#        """


rule compute_diamond_index:
    input: KEGG_DB
    output: DIAMOND_OUT+'.dmnd'
   
    conda:
        "envs/diamond.yaml"
    params:
        db = DIAMOND_OUT
    shell:
        """
        diamond makedb --in {input} --db {params.db}
        """

rule prodigal:
    input: 'genomes/{id}.fa'
    output:
        faa='proteins/{id}.faa',
        gff='proteins/{id}_prodigal.gff'
    conda: "envs/prodigal.yaml"
    shell: 'prodigal -p meta -a {output.faa} -q -i {input} -f gff -o {output.gff}'

#rule run_blastx:
#    input: 'genomes/{id}.fa'
#    output: 'output/blast/{id}.blastx.tab'
#    params:
#        db = BLASTPREFIX
#    conda:
#        'envs/blast.yaml'
#    shell:
#        """
#        blastx -db {params.db} -query {input} -out {output} -outfmt 6 -evalue 1e-5
#        """

#rule run_blastp: 
#    input: 'proteins/{id}.faa'
#    output: 'output/blast/{id}.blastp.tab'
#    params:
#        db = BLASTPREFIX
#    conda:
#        'envs/blast.yaml'
#    shell:
#        """
#        blastp -db {params.db} -query {input} -out {output} -outfmt 6 -evalue 1e-5
#        """


rule diamond_map:
    input:
        dmnd = DIAMOND_OUT + '.dmnd',
        pep = 'proteins/{id}.faa',
    output: 'output/blast/{id}.blastp.tab' 
    conda:
        "envs/diamond.yaml"
    params:
        other="--outfmt 6 --sensitive --block-size 2.0 -p 8"
    shell:
        """
        diamond blastp --db {input.dmnd} -q {input.pep} -o {output} {params.other} 
        """ 


rule kegg_annot:
    input: 'output/blast/{id}.blastp.tab'
    output: 'output/kegg_table/{id}.kegg_table.tsv'
    params: 
        kegg = KEGG_DIR
    conda:
        'envs/kegg-annot.yaml'
    log:
        'logs/{id}.kegg_table.log'
    shell:
        """
        keggannot_genes2ko {params.kegg} {input} -o {output}  > {log} 
        """
