import json
import pandas as pd
from pathlib import Path

infile = Path("metadata/ncbi_assembly_metadata.jsonl")
out = Path("metadata/assembly_biosample_metadata.tsv")

rows = []

def pick(d, *keys):
    """
    Return first available key from dictionary.
    """
    if not isinstance(d, dict):
        return ""
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return ""

def nested(d, path):
    """
    Safely access nested dictionary path.
    """
    cur = d
    for p in path:
        if isinstance(cur, dict):
            cur = cur.get(p, {})
        else:
            return ""
    return cur if cur not in [None, {}] else ""

def collect_all_attributes(obj):
    """
    Recursively collect all BioSample-style attributes from any part of JSON.
    """
    attrs = {}

    if isinstance(obj, dict):
        if "attributes" in obj and isinstance(obj["attributes"], list):
            for a in obj["attributes"]:
                if isinstance(a, dict):
                    name = str(a.get("name", "")).lower().replace(" ", "_").replace("-", "_")
                    value = a.get("value", "")
                    if name:
                        attrs[name] = value

        for v in obj.values():
            attrs.update(collect_all_attributes(v))

    elif isinstance(obj, list):
        for item in obj:
            attrs.update(collect_all_attributes(item))

    return attrs

with open(infile) as f:
    for line in f:
        if not line.strip():
            continue

        j = json.loads(line)

        assembly_info = pick(j, "assembly_info", "assemblyInfo")
        organism_block = pick(j, "organism")
        biosample = pick(assembly_info, "biosample", "bioSample")

        attrs = collect_all_attributes(j)

        rows.append({
            "assembly_accession": pick(j, "accession", "assembly_accession", "assemblyAccession"),
            "organism": pick(organism_block, "organism_name", "organismName", "name"),
            "strain": pick(
                pick(organism_block, "infraspecific_names", "infraspecificNames"),
                "strain"
            ),
            "biosample_accession": pick(biosample, "accession"),
            "bioproject_accession": pick(assembly_info, "bioproject_accession", "bioprojectAccession"),
            "assembly_level": pick(assembly_info, "assembly_level", "assemblyLevel"),
            "assembly_name": pick(assembly_info, "assembly_name", "assemblyName"),
            "submitter": pick(assembly_info, "submitter"),
            "isolation_source": pick(attrs, "isolation_source", "isolation_source_material", "source_material_id"),
            "host": pick(attrs, "host"),
            "geo_loc_name": pick(attrs, "geo_loc_name", "geographic_location"),
            "collection_date": pick(attrs, "collection_date"),
            "sample_type": pick(attrs, "sample_type"),
            "source_type": pick(attrs, "source_type"),
            "disease": pick(attrs, "disease"),
            "host_disease": pick(attrs, "host_disease"),
            "clinical_source": pick(attrs, "clinical_source"),
            "lat_lon": pick(attrs, "lat_lon"),
            "env_broad_scale": pick(attrs, "env_broad_scale"),
            "env_local_scale": pick(attrs, "env_local_scale"),
            "env_medium": pick(attrs, "env_medium"),
            "all_biosample_attributes": "; ".join([f"{k}={v}" for k, v in sorted(attrs.items())])
        })

df = pd.DataFrame(rows)
df.to_csv(out, sep="\t", index=False)

print("Saved:", out)
print("Rows:", len(df))
print("Columns:", len(df.columns))
print(df[["assembly_accession", "organism", "strain", "biosample_accession", "isolation_source", "host", "geo_loc_name"]].head())
