#!/usr/bin/env python3
"""Read-only raw-data discovery. No biological inference or filtering is performed.

Run from the study root: python scripts/discover_inputs.py --root .
All source_row numbers are one-based data-record numbers (header excluded).
JSON emitted to .yaml is valid YAML 1.2 and requires no PyYAML dependency.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

MISSING = {"", "NA", "N/A", "NaN", "nan", "null", "None", "-"}
DB = "DB_DETECTED"
INFERRED = "ANNOTATION_INFERRED"
MIXED = f"{DB}|{INFERRED}"


def write_table(rows, path):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(path, sep="\t", index=False)


def read_table(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def make_aliases(metadata):
    """Never truncate accession versions except via a unique metadata lookup."""
    candidates = defaultdict(set)
    provenance = defaultdict(list)
    for row, rec in enumerate(metadata.to_dict("records"), 1):
        canonical = rec["assembly_accession"]
        aliases = [canonical, re.sub(r"\.\d+$", "", canonical)]
        if isinstance(rec.get("genome_file"), str):
            aliases.append(re.sub(r"\.(?:fna|fa|fasta)(?:\.gz)?$", "", rec["genome_file"]))
        for alias in aliases:
            candidates[alias].add(canonical)
            provenance[alias].append(row)
    aliases = {k: next(iter(v)) for k, v in candidates.items() if len(v) == 1}
    ambiguous = {k: sorted(v) for k, v in candidates.items() if len(v) > 1}
    return aliases, ambiguous, provenance


def resolve_sample(value, aliases):
    if value in aliases:
        return aliases[value]
    stripped = re.sub(r"(?:\.mapping(?:\.potential)?\.ARG\.deeparg|\.tsv\.amrfinderplus|\.txt\.rgi)$", "", value)
    return aliases.get(stripped)


def table_profile(path, relative, sample, evidence, headerless=False):
    """Full-file scan: actual pandas dtypes, missing tokens, widths and record counts."""
    names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"] if headerless else None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     header=None if headerless else 0, names=names)
    cols = []
    for col in df:
        s = df[col]
        missing = s.isin(MISSING)
        v = s[~missing]
        numeric = pd.to_numeric(v, errors="coerce")
        dtype = "all_missing" if len(v) == 0 else "string"
        if len(v) and numeric.notna().all():
            dtype = "integer" if (numeric % 1 == 0).all() else "float"
        cols.append(dict(source_file=relative, column=col, dtype=dtype, n_rows=len(df),
                         n_missing=int(missing.sum()), n_empty=int((s == "").sum()),
                         n_dash=int((s == "-").sum()), n_literal_NA=int((s == "NA").sum()),
                         n_unique_nonmissing=int(v.nunique()), evidence_class=evidence))
    width_counts = Counter()
    with path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        if not headerless:
            next(reader, None)
        for rec in reader:
            if rec:
                width_counts[len(rec)] += 1
    rec = dict(source_file=relative, bytes=path.stat().st_size, sha256=sha256(path),
               n_rows=len(df), n_columns=len(df.columns), sample_id=sample or "",
               delimiter="TAB", header_present=not headerless,
               irregular_width_rows=sum(n for w, n in width_counts.items() if w != len(df.columns)),
               evidence_class=evidence)
    head = dict(source_file=relative, columns=list(df.columns), first_data_records=df.head(3).to_dict("records"))
    return df, rec, cols, head


def write_run_summary(root):
    root = Path(root)
    summary = json.loads((root / "discovery/summary.json").read_text())
    coverage = pd.read_csv(root / "discovery/sample_layer_coverage.tsv", sep="\t")
    manifest = range(summary["files_profiled"])
    errors = range(summary["parse_errors"])
    unresolved = range(summary["unresolved_ids"])
    target = range(summary["clinical_metadata_n"])
    metadata = range(summary["metadata_n"])
    n_bad_bgc = summary["antismash_rows_with_outside_genes"]
    bgc_interval_audit = pd.read_csv(root / "discovery/antismash_interval_audit.tsv", sep="\t")
    arg_categories = {tuple(k.split(" / ")): v for k, v in summary["arg_source_categories"].items()}
    omitted = range(summary["amp_unclustered_rows"])
    rep_audit = pd.read_csv(root / "discovery/amp_representative_count_audit.tsv", sep="\t").to_dict("records")
    n_cov_over = summary["arg_coverage_over_100"]
    run_summary = ["# Run summary", "", "Status: Step 0 discovery completed. Downstream analysis paused for source clarification, as explicitly requested.", "",
        f"Primary files inspected: {len(manifest):,}. Parse errors: {len(errors)}. Unresolved input identifiers: {len(unresolved)}.",
        f"Clinical metadata cohort: {len(target):,}. Metadata records outside this cohort: {len(metadata)-len(target):,}; listed explicitly in discovery/metadata_outside_cohort.tsv.", "",
        "| Layer | Clinical genomes with records | Clinical genomes without records |", "|---|---:|---:|"]
    for layer in ["amp", "arg", "bgc", "cazyme"]:
        present = int((coverage.in_clinical_metadata & coverage[layer]).sum())
        absent = int((coverage.in_clinical_metadata & ~coverage[layer]).sum())
        run_summary.append(f"| {layer.upper()} | {present:,} | {absent:,} |")
    run_summary += ["", "No clinical genomes were dropped. Analysis thresholds have not been applied; tool-support >0 checks are inventory counts only.", "",
        "## Analyses not yet run", "",
        "QC filtering, layer biology, statistical testing, cross-layer integration, figures, XLSX export, the full CLI, the 20-genome smoke test and Study_Report.docx remain to be implemented after the Step 0 clarification checkpoint.",
        "No tree or pangenome matrix is present under the specified input root. Dedicated lineage/ST, assembly length/contig count and measured AST fields are absent from the two supplied metadata schemas. Dependent analyses require appropriate input; no synthetic values will be substituted.", "",
        "## Source warnings", "",
        f"- antiSMASH: {n_bad_bgc:,}/{len(bgc_interval_audit):,} rows list GFF-matched genes outside the stated interval; raw region output or a corrected export is needed.",
        "- RGI: Strict/Perfect/Loose confidence is not retained; original calls or ARG-run parameters are needed.",
        f"- DeepARG: {arg_categories.get(('deeparg','DeepARG potential ARG'),0):,} potential rows must remain separate from {arg_categories.get(('deeparg','DeepARG ARG'),0):,} ARG rows.",
        f"- AMP: {len(omitted):,} candidates lack supplied cluster membership; {sum(r['difference'] != 0 for r in rep_audit):,} representative counts disagree with cluster-table counts.",
        f"- RGI coverage exceeds 100% in {n_cov_over:,} records; preserve and flag these values.",
        "- AMP/BGC software versions and DRAMP/dbCAN database releases remain unverified. Parameters for inactive layers in a CAZyme run do not prove other runs' settings.",
        "", "Reproduction and available audit outputs are documented in README.md. The complete discovery implementation is scripts/discover_inputs.py."]
    (root / "RUN_SUMMARY.md").write_text("\n".join(run_summary) + "\n")


def main(root):
    root = root.resolve()
    out = root / "discovery"
    out.mkdir(exist_ok=True)
    (root / "config").mkdir(exist_ok=True)
    (root / "logs").mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", handlers=[
        logging.FileHandler(root / "logs/discovery.log", mode="w"), logging.StreamHandler()])
    metadata_path = "metadata/final_genome_metadata_classified.tsv"
    metadata = read_table(root / metadata_path)
    aliases, ambiguous, provenance = make_aliases(metadata)
    target = sorted(metadata.loc[metadata.isolate_class == "Clinical", "assembly_accession"])
    alias_rows = [dict(raw_id=a, canonical_id=b, source_file=metadata_path,
                       source_rows=";".join(map(str, sorted(set(provenance[a])))), evidence_class=MIXED)
                  for a, b in sorted(aliases.items())]
    write_table(alias_rows, out / "sample_aliases.tsv")
    (out / "ambiguous_aliases.json").write_text(json.dumps(ambiguous, indent=2) + "\n")
    excluded = metadata[metadata.isolate_class != "Clinical"].copy()
    excluded["source_row"] = excluded.index + 1
    excluded["source_file"] = metadata_path
    excluded["exclusion_reason"] = "outside requested Clinical cohort"
    excluded["evidence_class"] = MIXED
    write_table(excluded, out / "metadata_outside_cohort.tsv")
    m1 = pd.read_csv(root / metadata_path, sep="\t", keep_default_na=False)
    m2 = pd.read_csv(root / "metadata/assembly_biosample_metadata.tsv", sep="\t", keep_default_na=False)
    if m1.assembly_accession.duplicated().any() or m2.assembly_accession.duplicated().any():
        raise ValueError("Duplicate metadata accession: explicit reconciliation required")
    m2_rows = {r.assembly_accession: i + 1 for i, r in m2.iterrows()}
    m2 = m2.set_index("assembly_accession")
    metadata_disagreements = []
    for idx, rec in m1.iterrows():
        sid = rec.assembly_accession
        if sid not in m2.index:
            continue
        for column in sorted(set(m1.columns) & set(m2.columns)):
            v1, v2 = str(rec[column]), str(m2.loc[sid, column])
            if v1 != v2:
                metadata_disagreements.append(dict(sample_id=sid, column=column,
                    classified_value=v1, biosample_value=v2,
                    both_nonmissing=v1 not in MISSING and v2 not in MISSING,
                    in_clinical_cohort=rec.isolate_class == "Clinical",
                    source_file=metadata_path, source_row=idx + 1,
                    comparison_file="metadata/assembly_biosample_metadata.tsv", comparison_row=m2_rows[sid], evidence_class=MIXED))
    write_table(metadata_disagreements, out / "metadata_disagreements.tsv")

    schemas = {
        "amp": {"path": "1350_AMP/Ampcombi_summary.tsv", "sample_id": "sample_id", "gene_id": "CDS_id",
                "sequence": "aa_sequence", "contig": "contig_id", "start": "CDS_start", "end": "CDS_end",
                "probabilities": {"ampir": "prob_ampir", "amplify": "prob_amplify", "macrel": "prob_macrel"},
                "database_id": "DRAMP_ID", "identity": "pident", "query_coverage": "qcov",
                "evalue": "evalue", "hydrophobicity": "hydrophobicity", "isoelectric_point": "isoelectric_point"},
        "amp_clusters": {"path": "1350_AMP/Ampcombi_summary_cluster.tsv", "sample_id": "sample_id",
                         "gene_id": "CDS_id", "cluster_id": "cluster_id", "sequence": "aa_sequence"},
        "amp_representatives": {"path": "1350_AMP/Ampcombi_summary_cluster_representative_seq.tsv",
                                "representative_id": "seq_headers", "cluster_id": "index", "member_count": "total_cluster_members"},
        "arg": {"path": "1350_ARG/reports/hamronization_summarize/hamronization_combined_report.tsv",
                "sample_id": "input_file_name", "tool": "analysis_software_name", "gene": "gene_symbol",
                "gene_name": "gene_name", "drug_class": "drug_class", "mechanism": "resistance_mechanism",
                "identity": "sequence_identity", "coverage": "coverage_percentage", "contig_or_protein": "input_sequence_id",
                "start": "input_gene_start", "end": "input_gene_stop", "variation": "genetic_variation_type",
                "amino_acid_mutation": "amino_acid_mutation", "tool_version": "analysis_software_version",
                "database": "reference_database_name", "database_version": "reference_database_version"},
        "bgc": {"path": "1350_BGC/combgc_complete_summary.tsv", "sample_id": "sample_id",
                "tool": "Prediction_tool", "feature": "Product_class", "probability": "BGC_probability",
                "complete": "BGC_complete", "contig": "contig_id", "start": "BGC_start", "end": "BGC_end",
                "length": "BGC_length", "cds_count": "CDS_count", "gene_ids": "BGC_region_contig_ids", "mibig": "MIBiG_ID"},
        "cazyme": {"glob": "1350_cazyme/cazyme/dbcan/cazyme_annotation/*/*_overview.tsv",
                   "sample_id_from": "parent_directory", "gene_id": "Gene ID", "support": "#ofTools",
                   "tools": {"HMMER": "dbCAN_hmm", "dbCAN_sub": "dbCAN_sub", "DIAMOND": "DIAMOND"},
                   "recommended": "Recommend Results", "substrate": "Substrate"},
        "cgc": {"glob": "1350_cazyme/cazyme/dbcan/cgc/*/*_cgc_standard_out.tsv",
                "sample_id_from": "parent_directory", "cluster_id": "CGC#", "gene_type": "Gene Type",
                "contig": "Contig ID", "gene_id": "Protein ID", "start": "Gene Start", "end": "Gene Stop"},
        "substrate": {"glob": "1350_cazyme/cazyme/dbcan/substrate/*/*_substrate_prediction.tsv",
                      "sample_id_from": "parent_directory", "cluster_id": "#cgcid", "pul_id": "PULID",
                      "pul_substrate": "dbCAN-PUL substrate", "subfamily_substrate": "dbCAN-sub substrate"},
        "metadata": {"path": metadata_path, "sample_id": "assembly_accession", "filename_alias": "genome_file",
                     "cohort": "isolate_class", "cohort_value": "Clinical", "source": "isolation_source",
                     "country": "geo_loc_name", "collection_date": "collection_date", "assembly_level": "assembly_level"}}
    mapping = {"schema_version": 1, "validated_against": "discovery/column_inventory.tsv",
               "sample_aliases": "discovery/sample_aliases.tsv", "sources": schemas,
               "coordinate_notes": {
                   "arg": "DeepARG start/end are protein alignment coordinates, not genomic intervals. Exclude from genomic proximity unless mapped through gene annotation.",
                   "bgc": "antiSMASH length=end-start+1; DeepBGC and GECCO length=end-start in supplied summary. Verify raw coordinate origin before proximity.",
                   "amp": "Coordinate origin not established by table headers; require producer documentation or raw annotation validation.",
                   "cgc_gff": "GFF genomic coordinates are one-based inclusive; annotated maximum end is not contig length."}}
    # Validate every declared data-column mapping against actual headers, including all per-genome schema variants below.
    for name, spec in schemas.items():
        p = root / spec["path"] if "path" in spec else sorted(root.glob(spec["glob"]))[0]
        actual = set(pd.read_csv(p, sep="\t", nrows=0).columns)
        for key, val in spec.items():
            if key in {"path", "glob", "sample_id_from", "cohort_value"}:
                continue
            for column in (val.values() if isinstance(val, dict) else [val]):
                if column not in actual:
                    raise ValueError(f"{name}: absent mapped column {column}")
    (root / "config/schema_map.yaml").write_text(json.dumps(mapping, indent=2) + "\n")

    paths = []
    for name in ["amp", "amp_clusters", "amp_representatives", "arg", "bgc", "metadata"]:
        paths.append((name, root / schemas[name]["path"], "", MIXED, False))
    paths.append(("metadata_secondary", root / "metadata/assembly_biosample_metadata.tsv", "", MIXED, False))
    dbcan_root = root / "1350_cazyme/cazyme/dbcan"
    for folder in ["cazyme_annotation", "cgc", "substrate"]:
        for p in sorted((dbcan_root / folder).glob("*/*")):
            if p.is_file() and p.suffix in {".tsv", ".out", ".tc", ".gff"}:
                family = p.name.replace(p.parent.name, "{sample}", 1)
                evidence = DB if family == "{sample}_diamond.out" else INFERRED
                if family == "{sample}_overview.tsv":
                    evidence = MIXED
                paths.append((family, p, resolve_sample(p.parent.name, aliases), evidence, p.suffix == ".gff"))
    manifest, columns, errors, frames = [], [], [], {}
    schema_variants = defaultdict(Counter)
    counts = defaultdict(int)
    samples = defaultdict(set)
    source_rows_by_tool = []
    cazyme_support = Counter()
    bgc_interval_audit = []
    overview_audit = []
    with (out / "input_heads.jsonl").open("w") as heads:
        for i, (family, p, sample, evidence, headerless) in enumerate(paths):
            relative = str(p.relative_to(root))
            try:
                df, rec, col, head = table_profile(p, relative, sample, evidence, headerless)
            except (pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as e:
                errors.append(dict(source_file=relative, error=str(e), evidence_class=evidence))
                continue
            rec["file_family"] = family
            manifest.append(rec)
            columns.extend(col)
            heads.write(json.dumps(head) + "\n")
            schema_variants[family][tuple(df.columns)] += 1
            counts[family] += len(df)
            if sample:
                samples[family].add(sample)
            if family in schemas or family == "metadata_secondary":
                frames[family] = df
            if family == "{sample}_overview.tsv":
                support = pd.to_numeric(df["#ofTools"], errors="coerce")
                cazyme_support.update(support.value_counts().to_dict())
                present = (~df[["dbCAN_hmm", "dbCAN_sub", "DIAMOND"]].isin(MISSING)).sum(axis=1)
                for idx in df.index[present != support]:
                    overview_audit.append(dict(source_file=relative, source_row=idx + 1,
                                              reported_support=support[idx], observed_support=present[idx], evidence_class=MIXED))
            if family == "{sample}_cgc.gff":
                bgc = frames["bgc"]
                selected = bgc[(bgc.Prediction_tool == "antiSMASH") & (bgc.canonical_id == sample)]
                if len(selected):
                    ids = df.attributes.str.extract(r"(?:^|;)protein_id=([^;]+)")[0]
                    genes = df.assign(protein_id=ids).set_index("protein_id")
                    for idx, row in selected.iterrows():
                        listed = row.BGC_region_contig_ids.split(";")
                        matched = genes.loc[genes.index.isin(listed)]
                        inside = ((matched.seqid == row.contig_id) &
                                  (pd.to_numeric(matched.start) >= float(row.BGC_start)) &
                                  (pd.to_numeric(matched.end) <= float(row.BGC_end)))
                        bgc_interval_audit.append(dict(source_file=schemas["bgc"]["path"], source_row=idx + 1,
                            sample_id=sample, contig_id=row.contig_id, bgc_start=row.BGC_start, bgc_end=row.BGC_end,
                            reported_cds_count=row.CDS_count, listed_gene_ids=len(listed), matched_gene_ids=len(matched),
                            matched_genes_outside_interval=int((~inside).sum()), validation_file=relative,
                            evidence_class=DB))
            if family == "bgc":
                df["canonical_id"] = df.sample_id.map(lambda x: resolve_sample(x, aliases))
            if i % 500 == 0:
                logging.info("Profiled %s/%s files", i + 1, len(paths))
    write_table(manifest, out / "file_inventory.tsv")
    write_table(columns, out / "column_inventory.tsv")
    write_table(errors if errors else [{"source_file": "", "error": "none", "evidence_class": MIXED}], out / "parse_errors.tsv")
    write_table(bgc_interval_audit, out / "antismash_interval_audit.tsv")
    write_table(overview_audit if overview_audit else [{"source_file": "", "issue": "no support discrepancies", "evidence_class": MIXED}], out / "dbcan_support_audit.tsv")

    presence = {s: {"sample_id": s, "in_clinical_metadata": True, "evidence_class": MIXED} for s in target}
    unresolved = []
    tool_rows = []
    for layer, sidcol, toolcol in [("amp", "sample_id", None), ("arg", "input_file_name", "analysis_software_name"), ("bgc", "sample_id", "Prediction_tool")]:
        df = frames[layer]
        df["canonical_id"] = df[sidcol].map(lambda x: resolve_sample(x, aliases))
        for idx, row in df[df.canonical_id.isna()].iterrows():
            unresolved.append(dict(layer=layer, raw_id=row[sidcol], source_file=schemas[layer]["path"], source_row=idx+1, evidence_class=MIXED))
        observed = set(df.canonical_id.dropna())
        for sid in observed | set(target):
            presence.setdefault(sid, {"sample_id": sid, "in_clinical_metadata": False, "evidence_class": MIXED})[layer] = sid in observed
        if toolcol:
            for tool, part in df.groupby(toolcol):
                cls = DB if tool in {"abricate", "amrfinderplus", "antiSMASH"} else INFERRED
                # RGI confidence absent: do not assert curated strict/perfect evidence.
                tool_rows.append(dict(layer=layer, tool=tool, n_rows=len(part), n_samples=part.canonical_id.nunique(),
                    source_file=schemas[layer]["path"], filter_applied=f"{toolcol} == {tool}", evidence_class=cls,
                    caveat="RGI confidence unavailable; evidence assignment provisional" if tool == "rgi" else ""))
        else:
            for tool, column in schemas["amp"]["probabilities"].items():
                part = df[pd.to_numeric(df[column]) > 0]
                tool_rows.append(dict(layer=layer, tool=tool, n_rows=len(part), n_samples=part.canonical_id.nunique(),
                    source_file=schemas[layer]["path"], filter_applied=f"{column} > 0 (reported support, not an analysis cutoff)", evidence_class=INFERRED))
    observed = samples["{sample}_overview.tsv"]
    for sid in observed | set(target):
        presence.setdefault(sid, {"sample_id": sid, "in_clinical_metadata": False, "evidence_class": MIXED})["cazyme"] = sid in observed
    for tool, column in schemas["cazyme"]["tools"].items():
        rows = [r for r in columns if r["column"] == column and r["source_file"].endswith("_overview.tsv")]
        tool_rows.append(dict(layer="cazyme", tool=tool, n_rows=sum(r["n_rows"] - r["n_missing"] for r in rows),
                             n_samples=sum(r["n_rows"] > r["n_missing"] for r in rows),
                             source_file=schemas["cazyme"]["glob"], filter_applied=f"{column} not a missing token", evidence_class=DB if tool=="DIAMOND" else INFERRED))
    coverage = pd.DataFrame(presence.values()).fillna(False).sort_values("sample_id")
    write_table(coverage, out / "sample_layer_coverage.tsv")
    write_table(tool_rows, out / "tool_contributions.tsv")
    write_table(unresolved if unresolved else [{"layer": "all", "raw_id": "", "status": "all identifiers resolved", "evidence_class": MIXED}], out / "unresolved_sample_ids.tsv")

    a, ac, ar, b = (frames[x] for x in ["amp", "amp_clusters", "arg", "bgc"])
    cluster_keys = set(zip(ac.sample_id, ac.CDS_id))
    omitted = a[[key not in cluster_keys for key in zip(a.sample_id, a.CDS_id)]].copy()
    omitted["source_row"] = omitted.index + 1
    omitted["source_file"] = schemas["amp"]["path"]
    omitted["evidence_class"] = INFERRED
    write_table(omitted[["sample_id", "canonical_id", "CDS_id", "aa_sequence", "DRAMP_ID", "source_row", "source_file", "evidence_class"]], out / "amp_rows_without_cluster.tsv")
    representative_counts = frames["amp_representatives"]
    cluster_sizes = ac.groupby("cluster_id").size()
    rep_audit = []
    for idx, rec in representative_counts.iterrows():
        reported = int(rec.total_cluster_members)
        observed = int(cluster_sizes.get(rec["index"], 0))
        rep_audit.append(dict(cluster_id=rec["index"], reported_members=reported, observed_cluster_rows=observed,
            difference=observed-reported, source_file=schemas["amp_representatives"]["path"], source_row=idx + 1,
            comparison_file=schemas["amp_clusters"]["path"], filter_applied=f"cluster_id == {rec['index']}", evidence_class=INFERRED))
    write_table(rep_audit, out / "amp_representative_count_audit.tsv")
    ar["source_category"] = "other"
    ar.loc[ar.input_file_name.str.endswith(".mapping.ARG.deeparg"), "source_category"] = "DeepARG ARG"
    ar.loc[ar.input_file_name.str.endswith(".mapping.potential.ARG.deeparg"), "source_category"] = "DeepARG potential ARG"
    arg_categories = ar.groupby(["analysis_software_name", "source_category"]).size().to_dict()
    n_cov_over = int((pd.to_numeric(ar.coverage_percentage, errors="coerce") > 100).sum())
    n_mutations = int((~ar.amino_acid_mutation.isin(MISSING)).sum())
    n_bad_bgc = sum(x["matched_genes_outside_interval"] > 0 for x in bgc_interval_audit)

    versions = []
    for keys, part in ar.groupby(["analysis_software_name", "analysis_software_version", "reference_database_name", "reference_database_version"]):
        versions.append(dict(layer="ARG", tool=keys[0], tool_version=keys[1], database=keys[2], database_version=keys[3],
            source_file=schemas["arg"]["path"], source_rows=";".join(str(x+1) for x in part.index),
            evidence_class=DB if keys[0] in {"abricate","amrfinderplus"} else INFERRED))
    version_file = root / "1350_cazyme/pipeline_info/nf_core_funcscan_software_mqc_versions.yml"
    for n, line in enumerate(version_file.read_text().splitlines(), 1):
        match = re.match(r"\s+(\S[^:]*):\s*(.+)$", line)
        if match:
            versions.append(dict(layer="CAZyme", tool=match[1], tool_version=match[2], database="", database_version="unreported",
                source_file=str(version_file.relative_to(root)), source_rows=f"physical_line:{n}", evidence_class=INFERRED))
    ham_version_path = "1350_ARG/reports/hamronization_summarize/versions.yml"
    for lineno, line in enumerate((root / ham_version_path).read_text().splitlines(), 1):
        match = re.match(r"\s+hamronization:\s*(.+)$", line)
        if match:
            versions.append(dict(layer="ARG", tool="hAMRonization", tool_version=match[1], database="", database_version="not applicable",
                source_file=ham_version_path, source_rows=f"physical_line:{lineno}", evidence_class=MIXED))
    write_table(versions, out / "tool_database_versions.tsv")
    param_rows = []
    provenance_files = []
    for p in sorted((root / "1350_cazyme/pipeline_info").glob("params_*.json")):
        params = json.loads(p.read_text())
        for k, v in params.items():
            param_rows.append(dict(source_file=str(p.relative_to(root)), parameter=k, value=json.dumps(v),
                                   run_cazyme_screening=params.get("run_cazyme_screening"),
                                   note="CAZyme run only; inactive AMP/ARG/BGC defaults do not prove settings for other runs", evidence_class=MIXED))
    write_table(param_rows, out / "pipeline_parameters.tsv")
    trace_rows = []
    for p in sorted((root / "1350_cazyme/pipeline_info").glob("execution_trace*")):
        trace = read_table(p)
        for (name, status), part in trace.assign(process=trace["name"].str.replace(r" \(.*$", "", regex=True)).groupby(["process", "status"]):
            trace_rows.append(dict(source_file=str(p.relative_to(root)), process=name, status=status, n_tasks=len(part), evidence_class=MIXED))
    write_table(trace_rows, out / "execution_statuses.tsv")
    for p in sorted(root.glob("1350_*/**/*")):
        if p.is_file() and ("pipeline_info" in p.parts or p.suffix == ".log" or p.name == "versions.yml"):
            provenance_files.append(dict(source_file=str(p.relative_to(root)), bytes=p.stat().st_size, sha256=sha256(p), evidence_class=MIXED))
    write_table(provenance_files, out / "provenance_manifest.tsv")

    variants = []
    for family, items in sorted(schema_variants.items()):
        for header, number in items.items():
            variants.append(dict(file_family=family, n_files=number, columns=json.dumps(header), evidence_class=MIXED))
    write_table(variants, out / "schema_variants.tsv")
    summary = dict(clinical_metadata_n=len(target), metadata_n=len(metadata), files_profiled=len(manifest),
                   parse_errors=len(errors), irregular_files=sum(r["irregular_width_rows"] > 0 for r in manifest),
                   unresolved_ids=len(unresolved), layer_samples={l: int(coverage[l].sum()) for l in ["amp", "arg", "bgc", "cazyme"]},
                   raw_rows={l: len(frames[l]) for l in ["amp", "amp_clusters", "amp_representatives", "arg", "bgc"]},
                   amp_cluster_count=ac.cluster_id.nunique(), amp_unclustered_rows=len(omitted),
                   amp_representative_count_disagreements=sum(r["difference"] != 0 for r in rep_audit),
                   arg_source_categories={" / ".join(k): v for k,v in arg_categories.items()},
                   arg_coverage_over_100=n_cov_over, explicit_amino_acid_mutation_rows=n_mutations,
                   antismash_rows_with_outside_genes=n_bad_bgc,
                   dbcan_support_counts={str(k): int(v) for k,v in cazyme_support.items()},
                   dbcan_support_discrepancies=len(overview_audit),
                   metadata_disagreements=len(metadata_disagreements),
                   metadata_nonmissing_disagreements=sum(r["both_nonmissing"] for r in metadata_disagreements),
                   file_family_rows=dict(counts))
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    md = ["# Step 0 data inventory", "", "Generated by `python scripts/discover_inputs.py --root .`. Raw inputs are never modified.",
          "All counts below are computed from the supplied files. This is a discovery audit, not the final biological analysis.",
          "", "## Cohort and identifier harmonisation", "",
          f"Metadata contain {len(metadata):,} assemblies; {len(target):,} have isolate_class Clinical. The other {len(metadata)-len(target):,} are retained in the alias audit but are outside the requested cohort.",
          "Three clinical assemblies use strain filenames that map explicitly to GCF accessions through genome_file. GECCO versionless accessions resolve only through unique metadata matches. Accession versions are never blindly stripped.",
          f"Unresolved identifiers: {len(unresolved)}. Ambiguous metadata aliases: {len(ambiguous)}.", "",
          "| Layer | Input rows | Harmonised genomes with records | Clinical genomes without records |", "|---|---:|---:|---:|"]
    for l in ["amp", "arg", "bgc", "cazyme"]:
        n = len(frames[l]) if l in frames else counts["{sample}_overview.tsv"]
        missing = int((coverage.in_clinical_metadata & ~coverage[l]).sum())
        md.append(f"| {l.upper()} | {n:,} | {int(coverage[l].sum()):,} | {missing:,} |")
    md += ["", "Per-sample membership and all missing/extra samples are explicit in discovery/sample_layer_coverage.tsv. A missing record is not assumed to be a biological absence.",
           "", "## Full input inspection", "",
           f"Scanned {len(manifest):,} primary input files completely, including every per-genome dbCAN annotation, CGC GFF/table and substrate table. Parse errors: {len(errors)}; files with irregular field counts: {summary['irregular_files']}.",
           "Every file has a SHA-256, record count and delimiter audit in discovery/file_inventory.tsv. Every actual column has an inferred dtype, missing-token counts and distinct-value count in discovery/column_inventory.tsv. The first three records of each file are in discovery/input_heads.jsonl. Headers are tab separated; GFF files are headerless. Missingness distinguishes blank, dash and literal NA; zero is never treated as missing.",
           "Discovery tables describe source provenance: an evidence_class value containing | explicitly lists the separate classes in that source; it never denotes consensus. Clinical metadata and software provenance are not biological evidence and inherit source scope only.",
           "", "| File family | Files | Rows | Schema variants |", "|---|---:|---:|---:|"]
    for fam in sorted(schema_variants):
        md.append(f"| {fam} | {sum(schema_variants[fam].values()):,} | {counts[fam]:,} | {len(schema_variants[fam])} |")
    md += ["", "## Observed schemas", "", "The machine-readable, header-validated mapping is config/schema_map.yaml (JSON syntax, valid YAML 1.2). Complete schemas follow; individual-file differences are in discovery/schema_variants.tsv."]
    for fam, variants_for_family in sorted(schema_variants.items()):
        md += ["", f"### {fam}"]
        for cols, n in variants_for_family.items():
            md += [f"{n:,} file(s):", "```text", "\t".join(cols), "```"]
    md += ["", "## Tools that actually contributed", "", "| Layer | Tool | Records with reported support | Genomes |", "|---|---|---:|---:|"]
    for row in tool_rows:
        md.append(f"| {row['layer']} | {row['tool']} | {row['n_rows']:,} | {row['n_samples']:,} |")
    md += ["", "fARGene is not represented. Do not infer execution from defaults in a different layer's parameter file.",
           "", "## Issues requiring resolution or explicit analysis policy", "",
           f"1. **antiSMASH region integrity:** {len(bgc_interval_audit):,} summary rows were cross-checked against matching protein IDs in dbCAN GFFs; {n_bad_bgc:,} rows list genes outside their stated BGC interval. The GFF check establishes a discrepancy, not its upstream cause. Per-row evidence is in discovery/antismash_interval_audit.tsv. Do not use these rows to estimate complete antiSMASH repertoires or known-cluster absence before inspecting original antiSMASH region JSON/GenBank output.",
           f"2. **DeepARG tiers:** {arg_categories.get(('deeparg','DeepARG potential ARG'),0):,} rows originate from mapping.potential.ARG and {arg_categories.get(('deeparg','DeepARG ARG'),0):,} from mapping.ARG. These strata must not be pooled as equivalent evidence. No DeepARG probability or coverage is retained in the harmonised report.",
           "3. **RGI confidence:** strict/perfect/loose confidence is absent from the harmonised columns. DB_DETECTED classification cannot be established for all RGI calls from this file alone. Original RGI files or ARG-run parameters are required to establish inclusion of loose hits.",
           f"4. **ARG mutations:** {n_mutations:,} records have explicit amino_acid_mutation text. Protein-variant labels alone are not proof of a reported specific mutation. AMRFinderPlus organism mode cannot be established from the supplied CAZyme-run parameters. Report this as unverified, not as definitely disabled.",
           f"5. **Coverage semantics:** {n_cov_over:,} RGI coverage values exceed 100%. Preserve and flag these rather than silently capping them. DeepARG start/end values are protein alignment coordinates and cannot be used directly for genomic proximity.",
           f"6. **AMP clustering:** {len(a):,} original candidates versus {len(ac):,} clustered rows, spanning {ac.cluster_id.nunique():,} reported clusters. {len(omitted):,} candidate rows have no supplied cluster assignment. Keep those rows unassigned; do not invent singleton families or infer exact-sequence clustering. The representative file contains identifiers/member counts, not amino-acid sequences. Sequences must be joined from the candidate table.",
           f"   Representative-file member counts disagree with grouped cluster-table row counts for {sum(r['difference'] != 0 for r in rep_audit):,} clusters; see discovery/amp_representative_count_audit.tsv. Recompute genome prevalence and row counts directly from cluster membership instead of trusting the representative total_cluster_members field.",
           f"7. **dbCAN support:** observed #ofTools counts are {dict(cazyme_support)}; {len(overview_audit):,} row(s) disagree with nonmissing tool columns. Two HMM-derived callers share related evidence; the requested support-based consensus label should not be described as statistical independence. Family-level agreement needs explicit parsing rather than equating gene-level support with family agreement.",
           "8. **Known BGCs:** MIBiG_ID is entirely missing. Broad NRPS/NRP classes do not establish pyoverdine, pyochelin or other specific products. antiSMASH PFAM lists are particularly unsafe for locus assignment until the interval discrepancy is resolved.",
           "9. **Proximity null:** annotation maximum coordinates are not contig lengths. A within-contig permutation null requires actual contig lengths and verified coordinate conventions.",
           "", "## Tool and database provenance", "",
           "CAZyme pipeline_info reports nf-core/funcscan v4.0.0-gaee3dc9, Nextflow 26.04.4, rundbcan 5.2.9 and Pyrodigal 3.6.3. These versions are scoped to that run. ARG reports ABRicate 1.0.1 / ncbi 2021-Mar-27; AMRFinderPlus 4.0.23 / database 2025-07-16.1; DeepARG 1.0.4 / database 2; RGI 6.0.5 / CARD 4.0.1. hAMRonization summarise reports 1.1.9.",
           "AMP predictor versions, AMPcombi/MMseqs versions, DRAMP release, BGC tool versions and dbCAN database release identifiers are not established by the provided version files. Do not fill these from current software releases or file timestamps. Local database files exist but their presence alone does not establish the database release.",
           "discovery/tool_database_versions.tsv links reported versions to exact source rows/lines; discovery/pipeline_parameters.tsv preserves every supplied parameter and its run-specific scope. provenance_manifest.tsv fingerprints supplied run/log files.",
           "", "## Metadata and optional analyses", "",
           "Both metadata tables are inspected and retained. The classified table establishes the clinical cohort. Country/source/year may support exploratory enrichment after normalising missing tokens and reporting denominators. All clinical assemblies are labelled Complete Genome in metadata; assembly length and contig-count covariates are not dedicated columns. No dedicated ST, lineage or measured AST columns are supplied. Free-text disease/host fields are not AST measurements.",
           f"The metadata tables differ at {len(metadata_disagreements):,} cells ({sum(r['both_nonmissing'] for r in metadata_disagreements):,} with both values outside the basic missing-token set). discovery/metadata_disagreements.tsv preserves both values and source rows; missing-like free text still needs normalisation. No conflicting value is silently overwritten.",
           "No phylogeny or Panaroo/Roary presence/absence matrix was found within the specified input root. Do not silently import other project folders or treat strain names as phylogenetic groups. Supervised ML remains optional and should not be enabled solely because country labels exist.",
           "", "## Analysis status", "",
           "Step 0 complete. Biological analysis and report generation have not started. The antiSMASH discrepancy and RGI confidence provenance require clarification before dependent scientific claims can be made. No results have been fabricated to fill these gaps."]
    (root / "00_data_inventory.md").write_text("\n".join(md) + "\n")
    write_run_summary(root)
    logging.info("DISCOVERY COMPLETE: %s", json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    main(args.root)
