#!/usr/bin/env python3
"""Build a machine-readable GEO/library manifest for GSE324375."""

from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
XML = ROOT / "metadata" / "GSE324375_family.xml"
OUT = ROOT / "metadata" / "geo_sample_manifest.tsv"
NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}


def clean(value: str | None) -> str:
    return " ".join((value or "").split())


tree = ET.parse(XML)
rows: list[dict[str, str]] = []

for sample in tree.findall("m:Sample", NS):
    title = sample.findtext("m:Title", default="", namespaces=NS)
    accession = sample.findtext("m:Accession", default="", namespaces=NS)
    channel = sample.find("m:Channel", NS)
    characteristics = {
        element.attrib.get("tag", ""): clean(element.text)
        for element in channel.findall("m:Characteristics", NS)
    }
    description = clean(sample.findtext("m:Description", default="", namespaces=NS))
    library_match = re.search(r"Library name:\s*(\S+)", description)
    supplementary = sample.find("m:Supplementary-Data", NS)
    relations = sample.findall("m:Relation", NS)
    relation_map = {r.attrib.get("type", ""): r.attrib.get("target", "") for r in relations}

    title_lower = title.lower()
    if title_lower.endswith(", rna"):
        modality = "GEX"
    elif title_lower.endswith(", adt"):
        modality = "ADT"
    elif title_lower.endswith(", hto"):
        modality = "HTO"
    else:
        modality = "unknown"

    rows.append(
        {
            "geo_accession": accession,
            "title": clean(title),
            "modality": modality,
            "library_name": library_match.group(1) if library_match else "",
            "source": clean(channel.findtext("m:Source", default="", namespaces=NS)),
            "tissue": characteristics.get("tissue", ""),
            "cell_line": characteristics.get("cell line", ""),
            "cell_type": characteristics.get("cell type", ""),
            "treatment_geo": characteristics.get("treatment", ""),
            "molecule": clean(channel.findtext("m:Molecule", default="", namespaces=NS)),
            "supplementary_type": supplementary.attrib.get("type", "") if supplementary is not None else "",
            "supplementary_url": clean(supplementary.text) if supplementary is not None else "",
            "biosample": relation_map.get("BioSample", "").rstrip("/").split("/")[-1],
            "sra_experiment": relation_map.get("SRA", "").split("term=")[-1],
        }
    )

fields = list(rows[0])
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} GEO sample records to {OUT}")
for modality in sorted({row['modality'] for row in rows}):
    print(f"{modality}: {sum(row['modality'] == modality for row in rows)}")
