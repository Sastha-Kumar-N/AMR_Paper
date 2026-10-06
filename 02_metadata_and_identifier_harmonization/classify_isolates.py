import pandas as pd
from pathlib import Path

genome_table = Path("accessions/genome_accession_table.tsv")
metadata_table = Path("metadata/assembly_biosample_metadata.tsv")

out_table = Path("metadata/final_genome_metadata_classified.tsv")
summary_table = Path("metadata/isolate_classification_summary.tsv")

genomes = pd.read_csv(genome_table, sep="\t").fillna("")
meta = pd.read_csv(metadata_table, sep="\t").fillna("")

df = genomes.merge(meta, on="assembly_accession", how="left").fillna("")

clinical_host_terms = [
    "homo sapiens", "human", "patient"
]

clinical_source_terms = [
    "clinical", "blood", "sputum", "urine", "urinary", "wound", "burn",
    "pus", "abscess", "catheter", "cystic fibrosis", "cf isolate",
    "respiratory", "lung", "bronchoalveolar", "cornea", "eye", "ear",
    "bone", "skin", "sepsis", "bacteremia", "bacteraemia", "infection",
    "infected", "stool", "feces", "faeces", "cancer", "nosocomial",
    "icu", "multidrug-resistant", "mdr"
]

environmental_terms = [
    "soil", "compost", "water", "river", "lake", "seawater", "sea water",
    "marine", "ocean", "wastewater", "waste water", "sewage", "effluent",
    "sediment", "sludge", "rhizosphere", "plant", "root", "leaf",
    "phyllosphere", "seed", "agricultural", "farm", "industrial",
    "oil", "petroleum", "air", "biofilm", "environment", "activated sludge",
    "drain", "sink", "tap water", "dental unit waterline"
]

def norm(x):
    return str(x).lower().strip()

def classify(row):
    host_text = norm(row.get("host", ""))
    
    source_text = " ".join([
        norm(row.get("isolation_source", "")),
        norm(row.get("sample_type", "")),
        norm(row.get("source_type", "")),
        norm(row.get("clinical_source", "")),
        norm(row.get("env_broad_scale", "")),
        norm(row.get("env_local_scale", "")),
        norm(row.get("env_medium", "")),
        norm(row.get("host_disease", "")),
        norm(row.get("disease", ""))
    ])

    all_text = source_text + " " + norm(row.get("all_biosample_attributes", ""))

    all_text = source_text + " " + norm(row.get("all_biosample_attributes", ""))

    host_is_clinical = any(t in host_text for t in clinical_host_terms) or any(t in all_text for t in clinical_host_terms)
    clinical_source_hit = any(t in source_text for t in clinical_source_terms)
    environmental_hit = any(t in source_text for t in environmental_terms)

    # Rule 1: direct human host = Clinical
    if host_is_clinical:
        return "Clinical"

    # Rule 2: clear environmental source and no human host = Environmental
    if environmental_hit:
        return "Environmental"

    # Rule 3: clear clinical material in source fields = Clinical
    if clinical_source_hit:
        return "Clinical"

    # Rule 4: fallback using full attributes
    if any(t in all_text for t in environmental_terms):
        return "Environmental"

    if any(t in all_text for t in clinical_source_terms):
        return "Clinical"

    return "Unknown_Other"

def evidence(row):
    text = " ".join([
        norm(row.get("host", "")),
        norm(row.get("isolation_source", "")),
        norm(row.get("sample_type", "")),
        norm(row.get("source_type", "")),
        norm(row.get("clinical_source", "")),
        norm(row.get("env_broad_scale", "")),
        norm(row.get("env_local_scale", "")),
        norm(row.get("env_medium", "")),
        norm(row.get("host_disease", "")),
        norm(row.get("disease", "")),
        norm(row.get("all_biosample_attributes", ""))
    ])

    hits = []

    for term in clinical_host_terms:
        if term in text:
            hits.append("clinical_host:" + term)

    for term in clinical_source_terms:
        if term in text:
            hits.append("clinical_source:" + term)

    for term in environmental_terms:
        if term in text:
            hits.append("environmental:" + term)

    return "; ".join(hits[:12])

df["isolate_class"] = df.apply(classify, axis=1)
df["classification_evidence"] = df.apply(evidence, axis=1)

df.to_csv(out_table, sep="\t", index=False)

summary = df["isolate_class"].value_counts().reset_index()
summary.columns = ["isolate_class", "count"]
summary.to_csv(summary_table, sep="\t", index=False)

print("Saved:", out_table)
print("Saved:", summary_table)
print()
print(summary)
