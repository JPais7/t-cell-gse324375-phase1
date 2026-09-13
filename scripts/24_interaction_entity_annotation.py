#!/usr/bin/env python3
"""Build an auditable mouse protein entity table from UniProtKB annotations."""
from pathlib import Path
import re
from io import StringIO
import requests
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/phase2"; CACHE=ROOT/"metadata/uniprot_mouse_entity_annotations.tsv"
SRC="https://rest.uniprot.org/uniprotkb/search"
FIELDS="accession,reviewed,gene_primary,protein_name,cc_function,cc_subcellular_location,go_f,go_c,ft_signal,ft_transmem"

def download():
 url=SRC; params={"query":"organism_id:10090 AND reviewed:true","format":"tsv","fields":FIELDS,"size":500}; chunks=[]
 while url:
  r=requests.get(url,params=params,timeout=120); r.raise_for_status(); chunks.append(r.text if not chunks else "\n".join(r.text.splitlines()[1:]))
  url=r.links.get("next",{}).get("url"); params=None
 CACHE.parent.mkdir(parents=True,exist_ok=True); CACHE.write_text("\n".join(c.rstrip("\n") for c in chunks)+"\n")

if not CACHE.exists(): download()
u=pd.read_csv(CACHE,sep="\t",dtype=str).fillna("")
u=u[u["Gene Names (primary)"].ne("")].copy()
u["reviewed_rank"]=u["Reviewed"].eq("reviewed").astype(int)
u=u.sort_values(["Gene Names (primary)","reviewed_rank"],ascending=[True,False]).drop_duplicates("Gene Names (primary)")
src=pd.read_csv(ROOT/"results/mechanisms/mechanisms_with_temporal_support.tsv.gz",sep="\t",usecols=["ligand","receptor"])
genes=sorted({p for col in ["ligand","receptor"] for value in src[col].dropna().astype(str) for p in value.split("_")})
u=u.set_index("Gene Names (primary)")
# Resolve symbols absent from the reviewed primary-name export (aliases and unreviewed-only entries).
for gene in [g for g in genes if g not in u.index]:
 r=requests.get(SRC,params={"query":f"gene_exact:{gene} AND organism_id:10090","format":"tsv","fields":FIELDS,"size":5},timeout=60); r.raise_for_status()
 q=pd.read_csv(StringIO(r.text),sep="\t",dtype=str).fillna("")
 if len(q):
  q["reviewed_rank"]=q["Reviewed"].eq("reviewed").astype(int); q=q.sort_values("reviewed_rank",ascending=False); q.index=[gene]+[f"{gene}__{i}" for i in range(1,len(q))]; u=pd.concat([u,q.iloc[:1]])

def yn(text,*terms): return any(t in text for t in terms)
rows=[]
for gene in genes:
 if gene in u.index:
  z=u.loc[gene]; text=" ".join(z.astype(str)).lower(); reviewed=z["Reviewed"]=="reviewed"
  protein=z["Protein names"].lower(); mf=(z["Gene Ontology (molecular function)"]+" "+z["Protein names"]).lower(); membrane=bool(z["Transmembrane"] or yn(z["Subcellular location [CC]"].lower(),"cell membrane","plasma membrane"))
  secreted=bool(z["Signal peptide"] or yn(text,"secreted","extracellular space","extracellular region"))
  receptor=yn(mf,"receptor activity","signaling receptor") or gene.startswith(("Tnfrsf","Ccr","Cxcr"))
  cytokine=yn(text,"cytokine activity","cytokine") or gene.startswith(("Il","Tnfsf"))
  chemokine=yn(text,"chemokine activity","chemokine") or gene.startswith(("Ccl","Cxcl"))
  protease=yn(mf,"peptidase activity","endopeptidase activity") or "metalloproteinase domain-containing" in protein
  enzyme=yn(mf,"catalytic activity","kinase activity","oxidoreductase","transferase activity","hydrolase activity","isomerase activity","lyase activity","decarboxylase activity","synthase activity","[ec ") or protease
  tf=yn(text,"dna-binding transcription factor activity","transcription factor")
  metabolic=yn(text,"glycolytic process","metabolic process","pyruvate","glucose metabolic","glucose-6-phosphate isomerase") and enzyme
  adhesion=yn(text,"cell adhesion","extracellular matrix","integrin binding","collagen") or gene.startswith(("Itga","Itgb","Col","Lam")) or gene in {"Thy1","Cd44","Hspg2"}
  ligand=(secreted and (cytokine or chemokine or yn(text,"growth factor activity","receptor ligand activity"))) or (membrane and gene.startswith("Tnfsf"))
  shedding=protease and membrane
  if tf: et="transcription factor"
  elif metabolic: et="intracellular enzyme/metabolic regulator"
  elif shedding: et="membrane protease/sheddase"
  elif membrane and ligand: et="membrane ligand"
  elif secreted and ligand: et="canonical soluble ligand"
  elif gene.startswith(("Itga","Itgb")) or receptor: et="receptor"
  elif adhesion: et="adhesion/contact molecule"
  elif enzyme: et="enzyme"
  elif membrane: et="membrane-associated protein"
  elif secreted: et="secreted/extracellular protein"
  else: et="unknown/unvalidated"
  rows.append({"gene":gene,"species":"Mus musculus","entity_type":et,"molecular_function":z["Gene Ontology (molecular function)"],"cellular_localization":z["Subcellular location [CC]"] or z["Gene Ontology (cellular component)"],"membrane_associated":membrane,"secreted":secreted,"ligand_supported":ligand,"receptor_supported":receptor,"enzyme":enzyme,"protease":protease,"transcription_factor":tf,"metabolic_regulator":metabolic,"cytokine":cytokine,"chemokine":chemokine,"adhesion_molecule":adhesion,"shedding_or_processing_factor":shedding,"evidence_source":f"UniProtKB {z['Entry']}; release retrieved by REST API","evidence_confidence":"HIGH" if reviewed else "MODERATE","biological_notes":z["Function [CC]"][:1000]})
 else:
  rows.append({"gene":gene,"species":"Mus musculus","entity_type":"unknown/unvalidated","molecular_function":"","cellular_localization":"","membrane_associated":False,"secreted":False,"ligand_supported":False,"receptor_supported":False,"enzyme":False,"protease":False,"transcription_factor":False,"metabolic_regulator":False,"cytokine":False,"chemokine":False,"adhesion_molecule":False,"shedding_or_processing_factor":False,"evidence_source":"No exact primary-gene match in cached UniProtKB mouse export","evidence_confidence":"LOW","biological_notes":""})
pd.DataFrame(rows).to_csv(OUT/"interaction_entity_annotation.tsv",sep="\t",index=False)
print(f"Annotated {len(rows)} mouse entities; high/moderate-confidence UniProt matches={sum(r['evidence_confidence']!='LOW' for r in rows)}")
