#!/usr/bin/env python3
"""Generate the conservative external perturbation evidence registry."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/phase2"
cols=["candidate","external_dataset","species","cell_type","context","perturbation_type","perturbed_gene","phenotype","effect","direction","FDR","dataset_match_quality","external_perturbation_support","external_context","source_url"]
rows=[
 ["Pkm","GSE289772","Mus musculus","CD8 T cell","activation and adoptive tumor therapy","pharmacologic PKM2 activation","Pkm","increased effector/recall and antitumor functions",None,"positive",None,"MODERATE_MATCH","CONTEXTUAL","Mouse CD8 context; pharmacologic activation, not gene-specific CRISPR and not validation of extracellular PKM-CD44","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289772"],
 ["Stat4","GSE314342","Homo sapiens","primary CD4 T cell","resting and stimulated","genome-scale CRISPRi","STAT4","dataset contains genome-scale perturbations; candidate-specific effect not extracted",None,None,None,"NOT_COMPARABLE","NOT_ASSESSED","Human CD4 context and no verified candidate-level effect in current integration","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342"],
 ["Adam10",None,None,None,None,None,"ADAM10","No directly comparable public T-cell/source-cell perturbation result verified",None,None,None,"NOT_COMPARABLE","NOT_ASSESSED","No eligible external perturbational evidence integrated",None],
]
pd.DataFrame(rows,columns=cols).to_csv(OUT/"external_perturbation_evidence.tsv",sep="\t",index=False)
