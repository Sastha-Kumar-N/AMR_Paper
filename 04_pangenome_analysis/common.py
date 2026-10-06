#!/usr/bin/env python3
"""Shared implementation for the whole-genome Panaroo publication analysis.

The functions in this module are intentionally conservative: original Panaroo
outputs are only read, all derived files are written under the requested output
directory, and expensive analyses have documented fallbacks.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import os
import platform
import random
import shutil
import statistics
import subprocess
import sys
import textwrap
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except Exception as exc:  # pragma: no cover
    raise RuntimeError("matplotlib is required for figure generation") from exc

try:
    import networkx as nx
except Exception as exc:  # pragma: no cover
    raise RuntimeError("networkx is required for network analysis") from exc

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

try:
    from scipy import sparse, stats
    from scipy.cluster import hierarchy
    from scipy.optimize import curve_fit
    from scipy.spatial.distance import squareform
    from scipy.sparse.linalg import svds
except Exception as exc:  # pragma: no cover
    raise RuntimeError("scipy is required for statistical analyses") from exc

try:
    from PIL import Image
except Exception:
    Image = None

try:
    from jinja2 import Template
except Exception:
    Template = None

try:
    from pyvis.network import Network
except Exception:
    Network = None


INPUT_FILES = [
    "combined_DNA_CDS.fasta",
    "combined_protein_CDS.fasta",
    "combined_protein_cdhit_out.txt",
    "combined_protein_cdhit_out.txt.clstr",
    "final_graph.gml",
    "pre_filt_graph.gml",
    "gene_data.csv",
    "gene_presence_absence.csv",
    "gene_presence_absence.Rtab",
    "gene_presence_absence_roary.csv",
    "struct_presence_absence.Rtab",
    "pan_genome_reference.fa",
    "summary_statistics.txt",
]

DIRS = [
    "00_audit",
    "01_scripts",
    "02_tables",
    "03_figures/individual",
    "03_figures/composite",
    "03_figures/figure_source_data",
    "04_networks/tables",
    "04_networks/graph_exports",
    "04_networks/cytoscape",
    "05_knowledge_graph/nodes",
    "05_knowledge_graph/edges",
    "05_knowledge_graph/graph_exports",
    "05_knowledge_graph/interactive",
    "06_interactive",
    "07_manuscript",
    "08_supplementary",
    "09_logs",
    "10_environment",
]

FREQ_ORDER = ["core", "soft_core", "shell", "cloud"]
FREQ_COLORS = {
    "core": "#0072B2",
    "soft_core": "#009E73",
    "shell": "#E69F00",
    "cloud": "#7F7F7F",
}
PALETTE = ["#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#D55E00", "#000000"]


@dataclass
class Context:
    input_dir: Path
    output_dir: Path
    threads: int = 20
    seed: int = 42
    permutations: int = 500

    def __post_init__(self) -> None:
        self.input_dir = self.input_dir.resolve()
        self.output_dir = self.output_dir.resolve()
        self.threads = max(1, min(int(self.threads), os.cpu_count() or 1, 20))
        random.seed(self.seed)
        np.random.seed(self.seed)


def build_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input-dir", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, default=Path("publication_deep_analysis"))
    parser.add_argument("--threads", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--cooccurrence-max-genes", type=int, default=1200)
    parser.add_argument("--cooccurrence-min-prevalence", type=int, default=20)
    parser.add_argument("--cooccurrence-max-prevalence-frac", type=float, default=0.95)
    parser.add_argument("--cooccurrence-min-shared", type=int, default=10)
    parser.add_argument("--accumulation-bootstrap", type=int, default=200)
    return parser


def context_from_args(args: argparse.Namespace) -> Context:
    return Context(args.input_dir, args.output_dir, args.threads, args.seed, args.permutations)


def ensure_dirs(ctx: Context) -> None:
    for rel in DIRS:
        (ctx.output_dir / rel).mkdir(parents=True, exist_ok=True)


def relpath(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except Exception:
        return str(path)


def log(ctx: Context, step: str, message: str) -> None:
    ensure_dirs(ctx)
    stamp = datetime.now().isoformat(timespec="seconds")
    line = f"{stamp}\t{step}\t{message}\n"
    with (ctx.output_dir / "09_logs/workflow.log").open("a", encoding="utf-8") as handle:
        handle.write(line)
    print(f"[{step}] {message}", flush=True)


def write_skip(ctx: Context, target: str, reason: str) -> None:
    path = ctx.output_dir / target
    path.parent.mkdir(parents=True, exist_ok=True)
    skip_path = path.with_suffix(path.suffix + ".skip_reason.txt") if path.suffix else Path(str(path) + ".skip_reason.txt")
    skip_path.write_text(reason.strip() + "\n", encoding="utf-8")
    log(ctx, "skip", f"{target}: {reason}")


def write_tsv(df: pd.DataFrame, path: Path, index: bool = False, **kwargs: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" or path.name.endswith(".tsv.gz") else None
    df.to_csv(path, sep="\t", index=index, compression=compression, **kwargs)


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def save_json(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str), encoding="utf-8")


def save_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def sha256_file(path: Path, block_size: int = 1024 * 1024 * 8) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detect_delimiter(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".rtab"} or path.name.endswith(".tsv.gz"):
        return "\\t"
    if suffix == ".csv":
        return ","
    try:
        with path.open("rb") as handle:
            sample = handle.read(8192).decode("utf-8", errors="replace")
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except Exception:
        return ""


def count_lines(path: Path) -> int | None:
    try:
        count = 0
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024 * 16), b""):
                count += chunk.count(b"\n")
        return count
    except Exception:
        return None


def get_header(path: Path, delimiter: str) -> list[str]:
    if not delimiter:
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t" if delimiter == "\\t" else delimiter)
            return next(reader)
    except Exception:
        return []


def fasta_count_headers(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                count += 1
    return count


def parse_summary_statistics(path: Path) -> dict[str, dict[str, Any]]:
    thresholds: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return thresholds
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            label, threshold, count = parts[0], parts[1], parts[2]
            key = label.lower().replace(" genes", "").replace(" ", "_")
            try:
                count_value = int(count)
            except ValueError:
                count_value = None
            thresholds[key] = {"label": label, "threshold_text": threshold, "count": count_value}
    return thresholds


def rtab_header(input_dir: Path, filename: str = "gene_presence_absence.Rtab") -> list[str]:
    with (input_dir / filename).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def load_rtab(ctx: Context, filename: str = "gene_presence_absence.Rtab") -> tuple[pd.Series, list[str], np.ndarray]:
    path = ctx.input_dir / filename
    if not path.exists():
        raise FileNotFoundError(path)
    log(ctx, "rtab", f"Loading {path.name}")
    df = pd.read_csv(path, sep="\t", dtype=str)
    genes = df.iloc[:, 0].astype(str)
    genomes = list(df.columns[1:])
    values = df.iloc[:, 1:].to_numpy(dtype=np.uint8, copy=True)
    return genes, genomes, values


def load_roary_annotations(ctx: Context) -> pd.DataFrame:
    path = ctx.input_dir / "gene_presence_absence_roary.csv"
    if not path.exists():
        path = ctx.input_dir / "gene_presence_absence.csv"
    if not path.exists():
        return pd.DataFrame()
    header = get_header(path, ",")
    metadata_cols = [c for c in header if not c.startswith(("GCA_", "GCF_"))]
    usecols = metadata_cols
    df = pd.read_csv(path, usecols=usecols, dtype=str)
    for col in ["No. isolates", "No. sequences", "Min group size nuc", "Max group size nuc", "Avg group size nuc"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def classify_prevalence(counts: np.ndarray | pd.Series, n_genomes: int) -> pd.Series:
    s = pd.Series(counts, copy=False)
    frac = s / float(n_genomes)
    out = pd.Series(np.full(len(s), "cloud", dtype=object), index=s.index)
    out[frac >= 0.15] = "shell"
    out[frac >= 0.95] = "soft_core"
    out[frac >= 0.99] = "core"
    return out


def is_uncharacterized_text(text: Any) -> bool:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return True
    t = str(text).strip().lower()
    if not t:
        return True
    markers = [
        "hypothetical",
        "uncharacterized",
        "unknown",
        "predicted protein",
        "putative protein",
        "conserved protein",
        "domain of unknown function",
    ]
    return any(m in t for m in markers)


def safe_label(text: Any, max_len: int = 60) -> str:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return ""
    t = " ".join(str(text).split())
    if len(t) <= max_len:
        return t
    return t[: max_len - 1] + "..."


def parse_gml_scalar(line: str) -> tuple[str, str] | None:
    stripped = line.strip()
    if not stripped or stripped in {"[", "]"}:
        return None
    if " " not in stripped:
        return None
    key, raw = stripped.split(" ", 1)
    raw = raw.strip()
    if raw.startswith('"') and raw.endswith('"'):
        raw = raw[1:-1]
    return key, raw


def cast_scalar(raw: str) -> Any:
    if raw == "":
        return ""
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw


def stream_gml_basic(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Parse a Panaroo GML file while discarding huge sequence attributes."""
    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    graph_attrs: dict[str, list[Any]] = defaultdict(list)
    mode: str | None = None
    current: dict[str, Any] = {}
    repeated: dict[str, int] = defaultdict(int)
    min_len: int | None = None
    max_len: int | None = None
    sum_len = 0
    n_len = 0
    skip_large = {"protein", "dna", "mergedDNA", "centroid", "seqIDs", "geneIDs", "members"}
    keep_node = {
        "id",
        "label",
        "size",
        "maxLenId",
        "hasEnd",
        "annotation",
        "description",
        "longCentroidID",
        "paralog",
        "genomeIDs",
        "degrees",
        "name",
    }
    keep_edge = {"source", "target", "size", "genomeIDs"}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            s = line.strip()
            if s == "node [":
                mode = "node"
                current = {}
                repeated = defaultdict(int)
                min_len = max_len = None
                sum_len = 0
                n_len = 0
                continue
            if s == "edge [":
                mode = "edge"
                current = {}
                repeated = defaultdict(int)
                continue
            if s == "]" and mode == "node":
                if n_len:
                    current["length_min"] = min_len
                    current["length_max"] = max_len
                    current["length_mean"] = sum_len / n_len
                    current["length_count"] = n_len
                current["members_count"] = repeated.get("members", 0)
                current["seqIDs_count"] = repeated.get("seqIDs", 0)
                current["geneIDs_count"] = repeated.get("geneIDs", 0)
                nodes.append(current)
                mode = None
                continue
            if s == "]" and mode == "edge":
                current["members_count"] = repeated.get("members", 0)
                edges.append(current)
                mode = None
                continue
            parsed = parse_gml_scalar(line)
            if parsed is None:
                continue
            key, raw = parsed
            if mode is None:
                if key == "isolateNames":
                    graph_attrs[key].append(raw)
                continue
            if key in skip_large:
                repeated[key] += 1
                continue
            if mode == "node" and key == "lengths":
                try:
                    val = int(raw)
                    min_len = val if min_len is None else min(min_len, val)
                    max_len = val if max_len is None else max(max_len, val)
                    sum_len += val
                    n_len += 1
                except ValueError:
                    pass
                continue
            if mode == "node" and key in keep_node:
                current[key] = cast_scalar(raw)
            elif mode == "edge" and key in keep_edge:
                current[key] = cast_scalar(raw)
    node_df = pd.DataFrame(nodes)
    edge_df = pd.DataFrame(edges)
    if "id" in node_df.columns:
        node_df["node_id"] = pd.to_numeric(node_df["id"], errors="coerce").astype("Int64")
    if "label" in node_df.columns:
        node_df["graph_label"] = node_df["label"].astype(str)
    if "name" in node_df.columns:
        node_df["gene_cluster_id"] = node_df["name"].astype(str)
    if "size" in node_df.columns:
        node_df["graph_prevalence"] = pd.to_numeric(node_df["size"], errors="coerce")
    if not edge_df.empty:
        for col in ["source", "target", "size"]:
            if col in edge_df.columns:
                edge_df[col] = pd.to_numeric(edge_df[col], errors="coerce")
    return node_df, edge_df, {k: v for k, v in graph_attrs.items()}


def graph_from_tables(node_df: pd.DataFrame, edge_df: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    for row in node_df.to_dict("records"):
        node_id = int(row.get("node_id") if not pd.isna(row.get("node_id")) else row.get("id"))
        attrs = {}
        for key in ["gene_cluster_id", "graph_label", "annotation", "description", "graph_prevalence", "length_mean", "length_min", "length_max"]:
            value = row.get(key)
            if value is not None and not (isinstance(value, float) and math.isnan(value)):
                attrs[key] = value
        graph.add_node(node_id, **attrs)
    for row in edge_df.to_dict("records"):
        try:
            u = int(row["source"])
            v = int(row["target"])
        except Exception:
            continue
        if u == v:
            continue
        weight = row.get("size", 1)
        if weight is None or pd.isna(weight):
            weight = 1
        graph.add_edge(u, v, weight=float(weight), size=float(weight))
    return graph


def parse_graph_cached(ctx: Context, graph_name: str = "final_graph.gml") -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    cache_nodes = ctx.output_dir / "04_networks/tables" / f"{Path(graph_name).stem}_parsed_nodes.tsv.gz"
    cache_edges = ctx.output_dir / "04_networks/tables" / f"{Path(graph_name).stem}_parsed_edges.tsv.gz"
    cache_attrs = ctx.output_dir / "04_networks/tables" / f"{Path(graph_name).stem}_graph_attrs.json"
    if cache_nodes.exists() and cache_edges.exists():
        return read_tsv(cache_nodes), read_tsv(cache_edges), json.loads(cache_attrs.read_text()) if cache_attrs.exists() else {}
    log(ctx, "gml", f"Streaming {graph_name}")
    node_df, edge_df, attrs = stream_gml_basic(ctx.input_dir / graph_name)
    write_tsv(node_df, cache_nodes)
    write_tsv(edge_df, cache_edges)
    save_json(attrs, cache_attrs)
    return node_df, edge_df, attrs


def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[valid] = tmp
    return out


def setup_plot() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(ctx: Context, fig: plt.Figure, basename: str, composite: bool = False) -> None:
    folder = ctx.output_dir / "03_figures" / ("composite" if composite else "individual")
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(folder / f"{basename}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def source_data(ctx: Context, name: str, df: pd.DataFrame) -> None:
    write_tsv(df, ctx.output_dir / "03_figures/figure_source_data" / name)


def copy_for_supplement(ctx: Context, src: Path, name: str | None = None) -> None:
    if src.exists():
        dst = ctx.output_dir / "08_supplementary" / (name or src.name)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def inventory_and_validate(ctx: Context) -> None:
    ensure_dirs(ctx)
    log(ctx, "audit", "Starting input inventory")
    rows: list[dict[str, Any]] = []
    for path in sorted(ctx.input_dir.iterdir()):
        if path.name == ctx.output_dir.name or path.name.startswith("."):
            continue
        if path.is_dir():
            for sub in sorted(path.rglob("*")):
                if sub.is_file():
                    rows.append(inventory_one_file(ctx, sub))
        elif path.is_file():
            rows.append(inventory_one_file(ctx, path))
    inventory = pd.DataFrame(rows)
    write_tsv(inventory, ctx.output_dir / "00_audit/input_file_inventory.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "00_audit/input_file_inventory.tsv", "Table_S1_input_file_inventory.tsv")

    log(ctx, "audit", "Recording Python environment")
    env_text = [
        f"analysis_timestamp\t{datetime.now(timezone.utc).isoformat()}",
        f"python_version\t{sys.version.replace(os.linesep, ' ')}",
        f"executable\t{sys.executable}",
        f"platform\t{platform.platform()}",
        f"uname\t{' '.join(platform.uname())}",
        f"cpu_count\t{os.cpu_count()}",
        f"threads_configured\t{ctx.threads}",
        f"seed\t{ctx.seed}",
        f"permutations\t{ctx.permutations}",
    ]
    if psutil:
        vm = psutil.virtual_memory()
        env_text.append(f"memory_total_bytes\t{vm.total}")
        env_text.append(f"memory_available_bytes\t{vm.available}")
    save_text("\n".join(env_text), ctx.output_dir / "10_environment/python_environment.txt")
    try:
        pip = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=freeze"],
            check=False,
            capture_output=True,
            text=True,
        )
        save_text(pip.stdout or pip.stderr, ctx.output_dir / "10_environment/requirements.txt")
    except Exception as exc:
        save_text(f"Unable to record pip freeze: {exc}", ctx.output_dir / "10_environment/requirements.txt")
    save_text(
        "name: panaroo_publication_deep_analysis\n"
        "channels:\n  - defaults\n"
        "dependencies:\n"
        f"  - python={sys.version_info.major}.{sys.version_info.minor}\n"
        "  - pip\n"
        "  - pip:\n"
        + "\n".join(f"      - {line}" for line in (ctx.output_dir / "10_environment/requirements.txt").read_text().splitlines() if line and not line.startswith("WARNING")),
        ctx.output_dir / "10_environment/environment.yml",
    )

    validation_rows: list[dict[str, Any]] = []
    consistency_rows: list[dict[str, Any]] = []
    unmatched_gene_rows: list[dict[str, Any]] = []
    unmatched_graph_rows: list[dict[str, Any]] = []
    warnings: list[str] = []

    summary = parse_summary_statistics(ctx.input_dir / "summary_statistics.txt")
    genes, genomes, matrix = load_rtab(ctx)
    n_genomes = len(genomes)
    n_genes = len(genes)
    prevalence = matrix.sum(axis=1)
    freq_class = classify_prevalence(prevalence, n_genomes)

    def add_check(check: str, severity: str, status: str, observed: Any = "", expected: Any = "", details: str = "") -> None:
        validation_rows.append(
            {
                "check": check,
                "severity": severity,
                "status": status,
                "observed": observed,
                "expected": expected,
                "details": details,
            }
        )

    add_check("rtab_genome_count", "info", "pass", n_genomes, "", "Genome columns in gene_presence_absence.Rtab")
    add_check("rtab_gene_cluster_count", "info", "pass", n_genes, summary.get("total", {}).get("count"), "Rows in gene_presence_absence.Rtab")
    dup_genes = int(genes.duplicated().sum())
    dup_genomes = len(genomes) - len(set(genomes))
    add_check("duplicate_gene_cluster_ids_rtab", "warning" if dup_genes else "info", "fail" if dup_genes else "pass", dup_genes)
    add_check("duplicate_genome_ids_rtab", "warning" if dup_genomes else "info", "fail" if dup_genomes else "pass", dup_genomes)
    nonbinary = int(((matrix != 0) & (matrix != 1)).sum())
    add_check("non_binary_rtab_values", "fatal" if nonbinary else "info", "fail" if nonbinary else "pass", nonbinary)
    empty_rows = int((prevalence == 0).sum())
    add_check("empty_gene_rows_rtab", "warning" if empty_rows else "info", "fail" if empty_rows else "pass", empty_rows)

    calc_counts = freq_class.value_counts().reindex(FREQ_ORDER, fill_value=0).to_dict()
    for cls in FREQ_ORDER:
        expected = summary.get(cls if cls != "soft_core" else "soft_core", {}).get("count")
        status = "pass" if expected is None or int(expected) == int(calc_counts[cls]) else "fail"
        sev = "warning" if status == "fail" else "info"
        add_check(f"{cls}_count_matches_summary_statistics", sev, status, calc_counts[cls], expected)

    rtab_set = set(genes)
    for table_name in ["gene_presence_absence.csv", "gene_presence_absence_roary.csv"]:
        path = ctx.input_dir / table_name
        if path.exists():
            ids = pd.read_csv(path, usecols=["Gene"], dtype=str)["Gene"]
            table_set = set(ids)
            consistency_rows.append({"source_a": "gene_presence_absence.Rtab", "source_b": table_name, "a_count": len(rtab_set), "b_count": len(table_set), "intersection": len(rtab_set & table_set), "a_not_b": len(rtab_set - table_set), "b_not_a": len(table_set - rtab_set)})
            for x in sorted(rtab_set - table_set)[:10000]:
                unmatched_gene_rows.append({"gene_cluster_id": x, "present_in": "gene_presence_absence.Rtab", "absent_from": table_name})
            for x in sorted(table_set - rtab_set)[:10000]:
                unmatched_gene_rows.append({"gene_cluster_id": x, "present_in": table_name, "absent_from": "gene_presence_absence.Rtab"})
            add_check(f"{table_name}_gene_cluster_count", "info", "pass", len(table_set), n_genes)

    struct_path = ctx.input_dir / "struct_presence_absence.Rtab"
    if struct_path.exists():
        struct_header = rtab_header(ctx.input_dir, "struct_presence_absence.Rtab")
        struct_rows = max(0, (count_lines(struct_path) or 1) - 1)
        add_check("struct_presence_absence_rows", "info", "pass", struct_rows, "")
        add_check("struct_presence_absence_genome_columns", "warning" if struct_header[1:] != genomes else "info", "pass" if struct_header[1:] == genomes else "fail", len(struct_header) - 1, n_genomes)

    for fasta_name in ["pan_genome_reference.fa"]:
        path = ctx.input_dir / fasta_name
        if path.exists():
            fasta_n = fasta_count_headers(path)
            add_check(f"{fasta_name}_sequence_count", "warning" if fasta_n != n_genes else "info", "pass" if fasta_n == n_genes else "fail", fasta_n, n_genes)

    node_df, edge_df, graph_attrs = parse_graph_cached(ctx, "final_graph.gml")
    pre_node_df, pre_edge_df, _ = parse_graph_cached(ctx, "pre_filt_graph.gml")
    add_check("final_graph_node_count", "warning" if len(node_df) != n_genes else "info", "pass" if len(node_df) == n_genes else "fail", len(node_df), n_genes)
    add_check("final_graph_edge_count", "info", "pass", len(edge_df))
    add_check("pre_filt_graph_node_count", "info", "pass", len(pre_node_df))
    add_check("pre_filt_graph_edge_count", "info", "pass", len(pre_edge_df))
    graph_names = set(node_df.get("gene_cluster_id", pd.Series(dtype=str)).dropna().astype(str))
    consistency_rows.append({"source_a": "gene_presence_absence.Rtab", "source_b": "final_graph.gml:name", "a_count": len(rtab_set), "b_count": len(graph_names), "intersection": len(rtab_set & graph_names), "a_not_b": len(rtab_set - graph_names), "b_not_a": len(graph_names - rtab_set)})
    for x in sorted(rtab_set - graph_names)[:10000]:
        unmatched_gene_rows.append({"gene_cluster_id": x, "present_in": "gene_presence_absence.Rtab", "absent_from": "final_graph.gml:name"})
    for x in sorted(graph_names - rtab_set)[:10000]:
        unmatched_graph_rows.append({"graph_node_gene_cluster_id": x, "present_in": "final_graph.gml:name", "absent_from": "gene_presence_absence.Rtab"})

    graph_attr_summary = []
    for col in node_df.columns:
        graph_attr_summary.append({"attribute_scope": "node", "attribute": col, "present_count": int(node_df[col].notna().sum())})
    for col in edge_df.columns:
        graph_attr_summary.append({"attribute_scope": "edge", "attribute": col, "present_count": int(edge_df[col].notna().sum())})
    write_tsv(pd.DataFrame(graph_attr_summary), ctx.output_dir / "00_audit/graph_attribute_inventory.tsv")

    gene_data_path = ctx.input_dir / "gene_data.csv"
    if gene_data_path.exists():
        gd_header = get_header(gene_data_path, ",")
        add_check("gene_data_schema", "info", "pass", ";".join(gd_header), "", "gene_data.csv is a per-CDS sequence table; no direct one-row-per-cluster field was assumed.")

    validation = pd.DataFrame(validation_rows)
    consistency = pd.DataFrame(consistency_rows)
    unmatched_gene = pd.DataFrame(unmatched_gene_rows)
    unmatched_graph = pd.DataFrame(unmatched_graph_rows)
    write_tsv(validation, ctx.output_dir / "00_audit/validation_report.tsv")
    write_tsv(consistency, ctx.output_dir / "00_audit/cross_file_consistency.tsv")
    write_tsv(unmatched_gene, ctx.output_dir / "00_audit/unmatched_gene_clusters.tsv")
    write_tsv(unmatched_graph, ctx.output_dir / "00_audit/unmatched_graph_nodes.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "00_audit/validation_report.tsv", "Table_S2_dataset_validation.tsv")
    audit_summary = {
        "genomes": n_genomes,
        "gene_clusters_rtab": n_genes,
        "frequency_class_counts": {k: int(v) for k, v in calc_counts.items()},
        "summary_statistics_thresholds": summary,
        "final_graph_nodes": int(len(node_df)),
        "final_graph_edges": int(len(edge_df)),
        "pre_filt_graph_nodes": int(len(pre_node_df)),
        "pre_filt_graph_edges": int(len(pre_edge_df)),
        "fatal_checks_failed": int(((validation["severity"] == "fatal") & (validation["status"] == "fail")).sum()),
        "warning_checks_failed": int(((validation["severity"] == "warning") & (validation["status"] == "fail")).sum()),
    }
    save_json(audit_summary, ctx.output_dir / "00_audit/audit_summary.json")
    save_text(
        f"""
        # Data Quality Report

        The analysis audited {len(inventory)} input files under `{ctx.input_dir}`.

        `gene_presence_absence.Rtab` contains {n_genes:,} gene-cluster rows and {n_genomes:,} genome columns. The Panaroo frequency thresholds were read from `summary_statistics.txt`: core `(99% <= strains <= 100%)`, soft core `(95% <= strains < 99%)`, shell `(15% <= strains < 95%)`, and cloud `(0% <= strains < 15%)`.

        Recalculated frequency-class counts were: core {calc_counts['core']:,}, soft core {calc_counts['soft_core']:,}, shell {calc_counts['shell']:,}, cloud {calc_counts['cloud']:,}.

        The final graph contains {len(node_df):,} parsed nodes and {len(edge_df):,} parsed edges. Graph-node names matched {len(rtab_set & graph_names):,} Rtab gene-cluster identifiers. Any non-matches are listed in `unmatched_gene_clusters.tsv` and `unmatched_graph_nodes.tsv`.

        `gene_data.csv` was treated as a per-CDS sequence table because its inspected schema is `{", ".join(get_header(gene_data_path, ",")) if gene_data_path.exists() else "not available"}`. It was not silently treated as a one-row-per-cluster table.
        """,
        ctx.output_dir / "00_audit/data_quality_report.md",
    )


def inventory_one_file(ctx: Context, path: Path) -> dict[str, Any]:
    rel = relpath(path, ctx.input_dir)
    log(ctx, "inventory", f"Checksumming {rel}")
    stat = path.stat()
    delim = detect_delimiter(path)
    header = get_header(path, delim)
    line_count = None
    row_count = None
    if path.suffix.lower() in {".csv", ".rtab", ".tsv", ".txt", ".gml", ".fa", ".fasta", ".clstr"} or path.name.endswith((".Rtab", ".fa")):
        line_count = count_lines(path)
        if line_count is not None and header:
            row_count = max(0, line_count - 1)
    return {
        "relative_path": rel,
        "size_bytes": stat.st_size,
        "modified_time_iso": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "sha256": sha256_file(path),
        "line_count": line_count,
        "row_count_if_delimited": row_count,
        "detected_delimiter": delim,
        "column_count": len(header),
        "columns": "|".join(header[:200]),
        "encoding_assumed": "utf-8 with replacement",
        "missing_value_conventions_observed": "empty string; NA not inferred during audit",
    }


def build_master_gene_table(ctx: Context) -> None:
    ensure_dirs(ctx)
    genes, genomes, matrix = load_rtab(ctx)
    prevalence = matrix.sum(axis=1).astype(int)
    n_genomes = len(genomes)
    master = pd.DataFrame(
        {
            "gene_cluster_id": genes.astype(str),
            "number_of_genomes_present": prevalence,
            "prevalence_fraction": prevalence / n_genomes,
            "prevalence_percentage": prevalence / n_genomes * 100,
        }
    )
    master["frequency_class"] = classify_prevalence(master["number_of_genomes_present"], n_genomes).values
    for cls in FREQ_ORDER:
        master[f"is_{cls}"] = master["frequency_class"].eq(cls)
    master["is_singleton"] = master["number_of_genomes_present"].eq(1)
    master["is_doubleton"] = master["number_of_genomes_present"].eq(2)

    roary = load_roary_annotations(ctx)
    if not roary.empty and "Gene" in roary.columns:
        roary = roary.rename(
            columns={
                "Gene": "gene_cluster_id",
                "Non-unique Gene name": "representative_gene_name",
                "Annotation": "annotation",
                "No. sequences": "number_of_gene_copies",
                "No. isolates": "roary_number_of_isolates",
                "Min group size nuc": "dna_length_min",
                "Max group size nuc": "dna_length_max",
                "Avg group size nuc": "dna_length_mean",
            }
        )
        keep = [c for c in roary.columns if c in {
            "gene_cluster_id",
            "representative_gene_name",
            "annotation",
            "number_of_gene_copies",
            "roary_number_of_isolates",
            "Genome Fragment",
            "Order within Fragment",
            "Accessory Fragment",
            "Accessory Order with Fragment",
            "QC",
            "dna_length_min",
            "dna_length_max",
            "dna_length_mean",
        }]
        master = master.merge(roary[keep], on="gene_cluster_id", how="left")

    node_df, edge_df, _ = parse_graph_cached(ctx, "final_graph.gml")
    graph = graph_from_tables(node_df, edge_df)
    deg = dict(graph.degree())
    wdeg = dict(graph.degree(weight="weight"))
    comp_map: dict[int, int] = {}
    comp_sizes: dict[int, int] = {}
    for idx, comp in enumerate(sorted(nx.connected_components(graph), key=len, reverse=True), start=1):
        comp_sizes[idx] = len(comp)
        for n in comp:
            comp_map[n] = idx
    try:
        core_numbers = nx.core_number(graph)
    except Exception:
        core_numbers = {n: 0 for n in graph.nodes}
    arts = set(nx.articulation_points(graph))
    bridges = list(nx.bridges(graph))
    bridge_nodes = {u for edge in bridges for u in edge}
    community_map = load_or_detect_communities(ctx, graph)
    node_metrics = pd.DataFrame(
        {
            "node_id": list(graph.nodes),
            "graph_degree": [deg[n] for n in graph.nodes],
            "weighted_degree": [wdeg[n] for n in graph.nodes],
            "connected_component": [comp_map.get(n, pd.NA) for n in graph.nodes],
            "connected_component_size": [comp_sizes.get(comp_map.get(n, 0), 0) for n in graph.nodes],
            "k_core_number": [core_numbers.get(n, 0) for n in graph.nodes],
            "articulation_point_status": [n in arts for n in graph.nodes],
            "bridge_related_status": [n in bridge_nodes for n in graph.nodes],
            "community_id": [community_map.get(n, pd.NA) for n in graph.nodes],
        }
    )
    graph_cols = node_df[["node_id", "gene_cluster_id", "graph_label", "annotation", "description", "graph_prevalence", "length_min", "length_max", "length_mean"]].copy()
    graph_cols = graph_cols.rename(columns={"annotation": "graph_annotation", "description": "graph_description", "length_mean": "graph_length_mean"})
    graph_table = graph_cols.merge(node_metrics, on="node_id", how="left")
    master = master.merge(graph_table, on="gene_cluster_id", how="left")
    master["product"] = master.get("annotation", pd.Series(index=master.index, dtype=object)).combine_first(master.get("graph_description", pd.Series(index=master.index, dtype=object)))
    master["functional_description"] = master["product"]
    master["sequence_length"] = master.get("dna_length_mean", pd.Series(index=master.index, dtype=float)).combine_first(master.get("graph_length_mean", pd.Series(index=master.index, dtype=float)))
    master["protein_length"] = np.where(master["sequence_length"].notna(), master["sequence_length"] / 3.0, np.nan)
    master["hypothetical_or_uncharacterized"] = master["product"].map(is_uncharacterized_text)
    master["available_graph_attributes"] = "final_graph.gml:name,label,size,annotation,description,lengths,degrees,genomeIDs"
    master = master.sort_values(["frequency_class", "number_of_genomes_present", "gene_cluster_id"], ascending=[True, False, True])
    out_tsv = ctx.output_dir / "02_tables/gene_cluster_master.tsv.gz"
    write_tsv(master, out_tsv)
    try:
        master.to_parquet(ctx.output_dir / "02_tables/gene_cluster_master.parquet", index=False)
    except Exception as exc:
        write_skip(ctx, "02_tables/gene_cluster_master.parquet", f"Parquet engine unavailable in current Python environment: {exc}")
    dd = data_dictionary(master, "gene_cluster_master")
    write_tsv(dd, ctx.output_dir / "02_tables/gene_cluster_master_data_dictionary.tsv")
    copy_for_supplement(ctx, out_tsv, "Table_S4_master_gene_cluster_annotations.tsv.gz")
    log(ctx, "master", f"Wrote {len(master):,} rows")


def data_dictionary(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    rows = []
    descriptions = {
        "gene_cluster_id": "Gene-cluster identifier from gene_presence_absence.Rtab / Panaroo table.",
        "number_of_genomes_present": "Number of genomes with presence value 1 in gene_presence_absence.Rtab.",
        "prevalence_fraction": "number_of_genomes_present divided by total genome count.",
        "frequency_class": "Core/soft_core/shell/cloud class using thresholds recorded in summary_statistics.txt.",
        "graph_degree": "Degree in final_graph.gml after self-loop removal.",
        "weighted_degree": "Sum of edge size attributes in final_graph.gml.",
        "community_id": "NetworkX Louvain community identifier when available.",
    }
    for col in df.columns:
        rows.append({"table": table_name, "column": col, "dtype": str(df[col].dtype), "description": descriptions.get(col, "Derived from supplied Panaroo files; see Methods for source and filtering.")})
    return pd.DataFrame(rows)


def load_master(ctx: Context) -> pd.DataFrame:
    path = ctx.output_dir / "02_tables/gene_cluster_master.tsv.gz"
    if not path.exists():
        build_master_gene_table(ctx)
    return read_tsv(path)


def pangenome_composition(ctx: Context) -> None:
    master = load_master(ctx)
    n = len(master)
    rows = []
    for cls in FREQ_ORDER:
        count = int((master["frequency_class"] == cls).sum())
        rows.append({"category": cls, "count": count, "percentage": count / n * 100})
    rows.append({"category": "core_plus_soft_core", "count": int(master["frequency_class"].isin(["core", "soft_core"]).sum()), "percentage": float(master["frequency_class"].isin(["core", "soft_core"]).mean() * 100)})
    rows.append({"category": "shell_plus_cloud", "count": int(master["frequency_class"].isin(["shell", "cloud"]).sum()), "percentage": float(master["frequency_class"].isin(["shell", "cloud"]).mean() * 100)})
    rows.append({"category": "singleton", "count": int(master["is_singleton"].sum()), "percentage": float(master["is_singleton"].mean() * 100)})
    rows.append({"category": "doubleton", "count": int(master["is_doubleton"].sum()), "percentage": float(master["is_doubleton"].mean() * 100)})
    comp = pd.DataFrame(rows)
    write_tsv(comp, ctx.output_dir / "02_tables/pangenome_composition.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/pangenome_composition.tsv", "Table_S3_whole_pangenome_composition.tsv")

    bins = [0, 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, master["number_of_genomes_present"].max()]
    labels = []
    for a, b in zip(bins[:-1], bins[1:]):
        labels.append(f"{a+1}-{b}" if a else "1")
    occ = master["number_of_genomes_present"].value_counts().sort_index().reset_index()
    occ.columns = ["genomes_present", "gene_cluster_count"]
    write_tsv(occ, ctx.output_dir / "02_tables/occupancy_distribution.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/occupancy_distribution.tsv", "Table_S6_gene_occupancy_distribution.tsv")
    binned = master.copy()
    binned["prevalence_bin"] = pd.cut(master["number_of_genomes_present"], bins=bins, labels=labels, include_lowest=True, duplicates="drop")
    gene_bins = binned.groupby("prevalence_bin", observed=True).size().reset_index(name="gene_cluster_count")
    write_tsv(gene_bins, ctx.output_dir / "02_tables/gene_frequency_bins.tsv")

    rare_defs = [
        ("singleton", master["number_of_genomes_present"].eq(1)),
        ("doubleton", master["number_of_genomes_present"].eq(2)),
        ("rare_le_5_genomes", master["number_of_genomes_present"].le(5)),
        ("rare_lt_1_percent", master["prevalence_fraction"].lt(0.01)),
        ("cloud_lt_15_percent", master["frequency_class"].eq("cloud")),
    ]
    rare = pd.DataFrame([{"definition": name, "count": int(mask.sum()), "percentage": float(mask.mean() * 100)} for name, mask in rare_defs])
    write_tsv(rare, ctx.output_dir / "02_tables/rare_gene_summary.tsv")

    ann_cov = (
        master.assign(annotated=~master["hypothetical_or_uncharacterized"].fillna(True))
        .groupby("frequency_class", observed=True)
        .agg(gene_clusters=("gene_cluster_id", "count"), annotated=("annotated", "sum"), uncharacterized=("hypothetical_or_uncharacterized", "sum"))
        .reset_index()
    )
    ann_cov["annotation_coverage_percentage"] = ann_cov["annotated"] / ann_cov["gene_clusters"] * 100
    write_tsv(ann_cov, ctx.output_dir / "02_tables/annotation_coverage_by_frequency_class.tsv")
    setup_plot()
    plot_composition_figures(ctx, comp, occ, rare, ann_cov, master)


def plot_composition_figures(ctx: Context, comp: pd.DataFrame, occ: pd.DataFrame, rare: pd.DataFrame, ann_cov: pd.DataFrame, master: pd.DataFrame) -> None:
    main = comp[comp["category"].isin(FREQ_ORDER)].copy()
    source_data(ctx, "pangenome_composition_source_data.tsv", main)
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.bar(main["category"], main["count"], color=[FREQ_COLORS[x] for x in main["category"]])
    ax.set_ylabel("Gene clusters")
    ax.set_title("Pangenome composition")
    ax.tick_params(axis="x", rotation=30)
    ax.yaxis.set_major_formatter(lambda x, pos: f"{int(x):,}")
    save_figure(ctx, fig, "pangenome_composition_bar")

    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    ax.pie(main["count"], labels=main["category"], colors=[FREQ_COLORS[x] for x in main["category"]], autopct="%1.1f%%", startangle=90, wedgeprops={"width": 0.45})
    ax.set_title("Frequency classes")
    save_figure(ctx, fig, "pangenome_composition_donut")

    source_data(ctx, "occupancy_distribution_source_data.tsv", occ)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(occ["genomes_present"], occ["gene_cluster_count"], color="#0072B2", lw=1)
    ax.set_xlabel("Genomes present")
    ax.set_ylabel("Gene clusters")
    ax.set_title("Gene-frequency spectrum")
    ax.yaxis.set_major_formatter(lambda x, pos: f"{int(x):,}")
    save_figure(ctx, fig, "gene_frequency_spectrum")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    low = occ[occ["genomes_present"] <= 100]
    ax.bar(low["genomes_present"], low["gene_cluster_count"], color="#7F7F7F")
    ax.set_yscale("log")
    ax.set_xlabel("Genomes present")
    ax.set_ylabel("Gene clusters (log scale)")
    ax.set_title("Low-frequency occupancy")
    save_figure(ctx, fig, "low_frequency_occupancy_log")

    sorted_prev = np.sort(master["number_of_genomes_present"].to_numpy())
    cdf = pd.DataFrame({"genomes_present": sorted_prev, "cumulative_fraction": np.arange(1, len(sorted_prev) + 1) / len(sorted_prev)})
    source_data(ctx, "cumulative_frequency_source_data.tsv", cdf)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(cdf["genomes_present"], cdf["cumulative_fraction"], color="#009E73")
    ax.set_xlabel("Genomes present")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Cumulative occupancy distribution")
    save_figure(ctx, fig, "cumulative_frequency_distribution")

    rank = pd.DataFrame({"rank": np.arange(1, len(sorted_prev) + 1), "genomes_present": sorted_prev[::-1]})
    source_data(ctx, "rank_prevalence_source_data.tsv", rank)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(rank["rank"], rank["genomes_present"], color="#D55E00", lw=1)
    ax.set_xlabel("Rank")
    ax.set_ylabel("Genomes present")
    ax.set_title("Rank-prevalence curve")
    save_figure(ctx, fig, "rank_prevalence_curve")

    source_data(ctx, "rare_gene_summary_source_data.tsv", rare)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.bar(rare["definition"], rare["count"], color="#CC79A7")
    ax.set_ylabel("Gene clusters")
    ax.set_title("Rare-cluster definitions")
    ax.tick_params(axis="x", rotation=35)
    ax.yaxis.set_major_formatter(lambda x, pos: f"{int(x):,}")
    save_figure(ctx, fig, "singleton_doubleton_rare_summary")

    source_data(ctx, "annotation_coverage_source_data.tsv", ann_cov)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    order = [x for x in FREQ_ORDER if x in set(ann_cov["frequency_class"])]
    vals = ann_cov.set_index("frequency_class").loc[order]
    ax.bar(vals.index, vals["annotation_coverage_percentage"], color=[FREQ_COLORS[x] for x in vals.index])
    ax.set_ylabel("Named/characterized clusters (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Annotation coverage by frequency class")
    save_figure(ctx, fig, "annotation_coverage_by_frequency_class")


def per_genome_gene_content(ctx: Context) -> None:
    master = load_master(ctx)
    genes, genomes, matrix = load_rtab(ctx)
    class_map = master.set_index("gene_cluster_id")["frequency_class"].to_dict()
    gene_classes = pd.Series(genes).map(class_map).fillna("unknown").to_numpy()
    singleton_mask = master.set_index("gene_cluster_id").reindex(genes)["is_singleton"].fillna(False).to_numpy(bool)
    rare_mask = master.set_index("gene_cluster_id").reindex(genes)["number_of_genomes_present"].fillna(0).le(5).to_numpy(bool)
    rows = []
    for j, genome in enumerate(genomes):
        col = matrix[:, j].astype(bool)
        total = int(col.sum())
        counts = {cls: int((col & (gene_classes == cls)).sum()) for cls in FREQ_ORDER}
        accessory = counts["shell"] + counts["cloud"]
        rows.append(
            {
                "genome_id": genome,
                "total_detected_gene_clusters": total,
                "core_gene_count": counts["core"],
                "soft_core_gene_count": counts["soft_core"],
                "shell_gene_count": counts["shell"],
                "cloud_gene_count": counts["cloud"],
                "accessory_gene_count": accessory,
                "rare_gene_count_le_5_genomes": int((col & rare_mask).sum()),
                "singleton_genes_unique_to_genome": int((col & singleton_mask).sum()),
                "percentage_accessory_content": accessory / total * 100 if total else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    med_total = df["total_detected_gene_clusters"].median()
    med_acc = df["accessory_gene_count"].median()
    mad_acc = stats.median_abs_deviation(df["accessory_gene_count"], scale="normal")
    df["deviation_from_dataset_median_total_genes"] = df["total_detected_gene_clusters"] - med_total
    df["deviation_from_dataset_median_accessory_genes"] = df["accessory_gene_count"] - med_acc
    df["standardized_accessory_gene_burden"] = (df["accessory_gene_count"] - med_acc) / mad_acc if mad_acc else 0.0
    q1, q3 = df["accessory_gene_count"].quantile([0.25, 0.75])
    iqr = q3 - q1
    df["iqr_accessory_outlier"] = (df["accessory_gene_count"] < q1 - 1.5 * iqr) | (df["accessory_gene_count"] > q3 + 1.5 * iqr)
    df["mad_accessory_outlier"] = df["standardized_accessory_gene_burden"].abs() > 3.5
    write_tsv(df, ctx.output_dir / "02_tables/per_genome_gene_content.tsv")
    out = df[df["iqr_accessory_outlier"] | df["mad_accessory_outlier"]].sort_values("standardized_accessory_gene_burden", key=lambda s: s.abs(), ascending=False)
    write_tsv(out, ctx.output_dir / "02_tables/per_genome_outliers.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/per_genome_gene_content.tsv", "Table_S5_per_genome_gene_content_statistics.tsv")
    setup_plot()
    plot_per_genome_figures(ctx, df)


def plot_per_genome_figures(ctx: Context, df: pd.DataFrame) -> None:
    source_data(ctx, "per_genome_gene_content_source_data.tsv", df)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.hist(df["total_detected_gene_clusters"], bins=40, color="#0072B2")
    ax.set_xlabel("Detected gene clusters per genome")
    ax.set_ylabel("Genomes")
    ax.set_title("Per-genome gene content")
    save_figure(ctx, fig, "distribution_total_genes_per_genome")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.hist(df["accessory_gene_count"], bins=40, color="#E69F00")
    ax.set_xlabel("Accessory clusters per genome")
    ax.set_ylabel("Genomes")
    ax.set_title("Accessory burden")
    save_figure(ctx, fig, "distribution_accessory_genes_per_genome")

    sorted_df = df.sort_values("accessory_gene_count").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    bottom = np.zeros(len(sorted_df))
    for cls, col in [("core", "core_gene_count"), ("soft_core", "soft_core_gene_count"), ("shell", "shell_gene_count"), ("cloud", "cloud_gene_count")]:
        ax.bar(np.arange(len(sorted_df)), sorted_df[col], bottom=bottom, color=FREQ_COLORS[cls], label=cls, linewidth=0)
        bottom += sorted_df[col].to_numpy()
    ax.set_xlabel("Genomes ranked by accessory burden")
    ax.set_ylabel("Gene clusters")
    ax.set_title("Frequency-class burden by genome")
    ax.legend(ncol=4, frameon=False)
    save_figure(ctx, fig, "stacked_per_genome_frequency_class_burden")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    rank = df.sort_values("total_detected_gene_clusters").reset_index(drop=True)
    ax.plot(np.arange(1, len(rank) + 1), rank["total_detected_gene_clusters"], color="#0072B2", lw=1)
    ax.set_xlabel("Genome rank")
    ax.set_ylabel("Detected gene clusters")
    ax.set_title("Ranked gene-content profile")
    save_figure(ctx, fig, "ranked_per_genome_gene_content")

    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.scatter(df["core_gene_count"] + df["soft_core_gene_count"], df["accessory_gene_count"], s=8, alpha=0.5, color="#009E73")
    ax.set_xlabel("Core + soft-core clusters")
    ax.set_ylabel("Accessory clusters")
    ax.set_title("Conserved vs accessory content")
    save_figure(ctx, fig, "core_versus_accessory_relationship")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.scatter(np.arange(len(df)), df["standardized_accessory_gene_burden"], s=8, alpha=0.5, color="#CC79A7")
    ax.axhline(3.5, color="black", lw=0.8, ls="--")
    ax.axhline(-3.5, color="black", lw=0.8, ls="--")
    ax.set_xlabel("Genome index")
    ax.set_ylabel("Robust z score")
    ax.set_title("Accessory-content outliers")
    save_figure(ctx, fig, "genome_level_outlier_plot")


def genome_bitsets(matrix: np.ndarray) -> list[int]:
    n_genes, n_genomes = matrix.shape
    bitsets: list[int] = []
    for j in range(n_genomes):
        bits = 0
        idx = np.flatnonzero(matrix[:, j])
        for i in idx:
            bits |= 1 << int(i)
        bitsets.append(bits)
    return bitsets


def accumulation_and_heaps(ctx: Context, bootstrap: int = 200) -> None:
    genes, genomes, matrix = load_rtab(ctx)
    n_genes, n_genomes = matrix.shape
    bitsets = genome_bitsets(matrix)
    rng = np.random.default_rng(ctx.seed)
    raw_path = ctx.output_dir / "02_tables/pangenome_accumulation_permutations.tsv.gz"
    log(ctx, "accumulation", f"Running {ctx.permutations} genome-order permutations with Python integer bitsets")
    with gzip.open(raw_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["permutation", "n_genomes", "pangenome_size", "core_genome_size", "new_gene_clusters", "singleton_clusters_in_cumulative_subset"])
        for p in range(ctx.permutations):
            order = rng.permutation(n_genomes)
            union = 0
            core = (1 << n_genes) - 1
            seen_once = 0
            seen_more = 0
            prev_pan = 0
            for step, genome_idx in enumerate(order, start=1):
                bits = bitsets[int(genome_idx)]
                union |= bits
                core &= bits
                newly_repeat = seen_once & bits
                seen_more |= newly_repeat
                seen_once = (seen_once | bits) & ~seen_more
                pan = union.bit_count()
                writer.writerow([p + 1, step, pan, core.bit_count(), pan - prev_pan, seen_once.bit_count()])
                prev_pan = pan
    raw = pd.read_csv(raw_path, sep="\t")
    summary = raw.groupby("n_genomes").agg(
        pangenome_mean=("pangenome_size", "mean"),
        pangenome_lower95=("pangenome_size", lambda x: np.quantile(x, 0.025)),
        pangenome_upper95=("pangenome_size", lambda x: np.quantile(x, 0.975)),
        core_mean=("core_genome_size", "mean"),
        core_lower95=("core_genome_size", lambda x: np.quantile(x, 0.025)),
        core_upper95=("core_genome_size", lambda x: np.quantile(x, 0.975)),
        new_gene_mean=("new_gene_clusters", "mean"),
        new_gene_lower95=("new_gene_clusters", lambda x: np.quantile(x, 0.025)),
        new_gene_upper95=("new_gene_clusters", lambda x: np.quantile(x, 0.975)),
        singleton_mean=("singleton_clusters_in_cumulative_subset", "mean"),
        singleton_lower95=("singleton_clusters_in_cumulative_subset", lambda x: np.quantile(x, 0.025)),
        singleton_upper95=("singleton_clusters_in_cumulative_subset", lambda x: np.quantile(x, 0.975)),
    ).reset_index()
    write_tsv(summary[["n_genomes", "pangenome_mean", "pangenome_lower95", "pangenome_upper95"]], ctx.output_dir / "02_tables/pangenome_accumulation_curve.tsv")
    write_tsv(summary[["n_genomes", "core_mean", "core_lower95", "core_upper95"]], ctx.output_dir / "02_tables/core_genome_decay_curve.tsv")
    write_tsv(summary[["n_genomes", "new_gene_mean", "new_gene_lower95", "new_gene_upper95"]], ctx.output_dir / "02_tables/new_gene_discovery_curve.tsv")
    write_tsv(summary, ctx.output_dir / "02_tables/singleton_accumulation_curve.tsv")

    params, boot = fit_heaps_models(raw, summary, ctx.seed, bootstrap)
    write_tsv(params, ctx.output_dir / "02_tables/heaps_law_model_parameters.tsv")
    write_tsv(boot, ctx.output_dir / "02_tables/heaps_law_bootstrap.tsv.gz")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/heaps_law_model_parameters.tsv", "Table_S7_heaps_law_model_results.tsv")
    setup_plot()
    plot_accumulation_figures(ctx, summary, params)


def power_law(n: np.ndarray, kappa: float, gamma: float) -> np.ndarray:
    return kappa * np.power(n, gamma)


def decay_law(n: np.ndarray, kappa: float, alpha: float) -> np.ndarray:
    return kappa * np.power(n, -alpha)


def model_stats(y: np.ndarray, yhat: np.ndarray, p: int) -> dict[str, float]:
    resid = y - yhat
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    n = len(y)
    rmse = math.sqrt(rss / n)
    r2 = 1 - rss / tss if tss else np.nan
    aic = n * math.log(rss / n) + 2 * p if rss > 0 else np.nan
    return {"rmse": rmse, "r_squared": r2, "aic": aic}


def fit_heaps_models(raw: pd.DataFrame, summary: pd.DataFrame, seed: int, bootstrap: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = summary["n_genomes"].to_numpy(float)
    y_pan = summary["pangenome_mean"].to_numpy(float)
    pan_params, pan_cov = curve_fit(power_law, x, y_pan, p0=(5000, 0.3), bounds=(0, np.inf), maxfev=10000)
    pan_pred = power_law(x, *pan_params)
    pan_stats = model_stats(y_pan, pan_pred, 2)
    x_new = summary.loc[summary["n_genomes"] >= 2, "n_genomes"].to_numpy(float)
    y_new = summary.loc[summary["n_genomes"] >= 2, "new_gene_mean"].to_numpy(float)
    new_params, new_cov = curve_fit(decay_law, x_new, y_new, p0=(1000, 0.5), bounds=(0, np.inf), maxfev=10000)
    new_pred = decay_law(x_new, *new_params)
    new_stats = model_stats(y_new, new_pred, 2)
    rows = [
        {"model": "pangenome_power_law_P(n)=kappa*n^gamma", "parameter": "kappa", "estimate": pan_params[0], "standard_error": math.sqrt(pan_cov[0, 0]) if np.isfinite(pan_cov[0, 0]) else np.nan, **pan_stats, "n_points": len(x), "permutations": raw["permutation"].nunique(), "interpretation": "Positive gamma indicates continuing pangenome expansion under this parameterization."},
        {"model": "pangenome_power_law_P(n)=kappa*n^gamma", "parameter": "gamma", "estimate": pan_params[1], "standard_error": math.sqrt(pan_cov[1, 1]) if np.isfinite(pan_cov[1, 1]) else np.nan, **pan_stats, "n_points": len(x), "permutations": raw["permutation"].nunique(), "interpretation": "Positive gamma indicates continuing pangenome expansion under this parameterization."},
        {"model": "new_gene_decay_Nnew(n)=kappa*n^-alpha", "parameter": "kappa", "estimate": new_params[0], "standard_error": math.sqrt(new_cov[0, 0]) if np.isfinite(new_cov[0, 0]) else np.nan, **new_stats, "n_points": len(x_new), "permutations": raw["permutation"].nunique(), "interpretation": "Larger alpha indicates faster decay in newly discovered genes per added genome."},
        {"model": "new_gene_decay_Nnew(n)=kappa*n^-alpha", "parameter": "alpha", "estimate": new_params[1], "standard_error": math.sqrt(new_cov[1, 1]) if np.isfinite(new_cov[1, 1]) else np.nan, **new_stats, "n_points": len(x_new), "permutations": raw["permutation"].nunique(), "interpretation": "Larger alpha indicates faster decay in newly discovered genes per added genome."},
    ]
    rng = np.random.default_rng(seed)
    boot_rows = []
    perms = raw["permutation"].unique()
    grouped = {p: g for p, g in raw.groupby("permutation")}
    for b in range(bootstrap):
        sample = rng.choice(perms, size=len(perms), replace=True)
        merged = pd.concat([grouped[int(p)] for p in sample], ignore_index=True)
        summ = merged.groupby("n_genomes").agg(pangenome_mean=("pangenome_size", "mean"), new_gene_mean=("new_gene_clusters", "mean")).reset_index()
        try:
            bp, _ = curve_fit(power_law, summ["n_genomes"].to_numpy(float), summ["pangenome_mean"].to_numpy(float), p0=pan_params, bounds=(0, np.inf), maxfev=10000)
            bn, _ = curve_fit(decay_law, summ.loc[summ["n_genomes"] >= 2, "n_genomes"].to_numpy(float), summ.loc[summ["n_genomes"] >= 2, "new_gene_mean"].to_numpy(float), p0=new_params, bounds=(0, np.inf), maxfev=10000)
            boot_rows.extend(
                [
                    {"bootstrap": b + 1, "model": "pangenome_power_law", "parameter": "kappa", "estimate": bp[0]},
                    {"bootstrap": b + 1, "model": "pangenome_power_law", "parameter": "gamma", "estimate": bp[1]},
                    {"bootstrap": b + 1, "model": "new_gene_decay", "parameter": "kappa", "estimate": bn[0]},
                    {"bootstrap": b + 1, "model": "new_gene_decay", "parameter": "alpha", "estimate": bn[1]},
                ]
            )
        except Exception:
            continue
    boot = pd.DataFrame(boot_rows)
    if not boot.empty:
        ci = boot.groupby(["model", "parameter"])["estimate"].quantile([0.025, 0.975]).unstack().reset_index().rename(columns={0.025: "bootstrap_ci95_low", 0.975: "bootstrap_ci95_high"})
        params = pd.DataFrame(rows)
        params["model_short"] = np.where(params["model"].str.startswith("pangenome_power"), "pangenome_power_law", "new_gene_decay")
        params = params.merge(ci, left_on=["model_short", "parameter"], right_on=["model", "parameter"], how="left", suffixes=("", "_ci"))
        params = params.drop(columns=[c for c in ["model_ci"] if c in params.columns])
    else:
        params = pd.DataFrame(rows)
    residuals = pd.DataFrame(
        {
            "n_genomes": np.concatenate([x, x_new]),
            "model": ["pangenome_power_law"] * len(x) + ["new_gene_decay"] * len(x_new),
            "observed": np.concatenate([y_pan, y_new]),
            "fitted": np.concatenate([pan_pred, new_pred]),
        }
    )
    residuals["residual"] = residuals["observed"] - residuals["fitted"]
    # The caller writes residual diagnostics source through the plotting function.
    params.attrs["residuals"] = residuals
    return params, boot


def plot_accumulation_figures(ctx: Context, summary: pd.DataFrame, params: pd.DataFrame) -> None:
    source_data(ctx, "accumulation_models_source_data.tsv", summary)
    x = summary["n_genomes"].to_numpy()
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(x, summary["pangenome_mean"], color="#0072B2")
    ax.fill_between(x, summary["pangenome_lower95"], summary["pangenome_upper95"], color="#0072B2", alpha=0.2)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Pangenome size")
    ax.set_title("Pangenome accumulation")
    save_figure(ctx, fig, "pangenome_accumulation_curve")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(x, summary["core_mean"], color="#009E73")
    ax.fill_between(x, summary["core_lower95"], summary["core_upper95"], color="#009E73", alpha=0.2)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Core-genome size")
    ax.set_title("Core-genome decay")
    save_figure(ctx, fig, "core_genome_decay_curve")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(x, summary["new_gene_mean"], color="#D55E00")
    ax.fill_between(x, summary["new_gene_lower95"], summary["new_gene_upper95"], color="#D55E00", alpha=0.2)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("New clusters per genome")
    ax.set_title("New-gene discovery")
    save_figure(ctx, fig, "new_gene_discovery_curve")

    pan = params[(params["model"].str.startswith("pangenome_power"))].set_index("parameter")["estimate"]
    fit = power_law(x.astype(float), pan["kappa"], pan["gamma"])
    fit_df = pd.DataFrame({"n_genomes": x, "observed_mean": summary["pangenome_mean"], "fitted": fit, "residual": summary["pangenome_mean"] - fit})
    source_data(ctx, "heaps_fit_source_data.tsv", fit_df)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(x, summary["pangenome_mean"], color="#0072B2", label="Observed mean")
    ax.plot(x, fit, color="black", ls="--", label="Power-law fit")
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Pangenome size")
    ax.set_title("Heaps' law fit")
    ax.legend(frameon=False)
    save_figure(ctx, fig, "heaps_law_observed_vs_fitted")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.scatter(fit_df["fitted"], fit_df["residual"], s=8, alpha=0.6, color="#CC79A7")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Fitted pangenome size")
    ax.set_ylabel("Residual")
    ax.set_title("Heaps' law residuals")
    save_figure(ctx, fig, "heaps_law_residual_diagnostic")

    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.4))
    axs[0].plot(x, summary["pangenome_mean"], color="#0072B2")
    axs[0].fill_between(x, summary["pangenome_lower95"], summary["pangenome_upper95"], color="#0072B2", alpha=0.2)
    axs[0].set_title("A. Pangenome")
    axs[1].plot(x, summary["core_mean"], color="#009E73")
    axs[1].fill_between(x, summary["core_lower95"], summary["core_upper95"], color="#009E73", alpha=0.2)
    axs[1].set_title("B. Core")
    axs[2].plot(x, summary["new_gene_mean"], color="#D55E00")
    axs[2].set_title("C. New clusters")
    for ax in axs:
        ax.set_xlabel("Genomes")
    save_figure(ctx, fig, "combined_accumulation_model_figure")


def annotation_landscape(ctx: Context) -> None:
    master = load_master(ctx)
    named = master[~master["hypothetical_or_uncharacterized"].fillna(True)].copy()
    top_named = named.sort_values("number_of_genomes_present", ascending=False).head(200)
    write_tsv(top_named[["gene_cluster_id", "representative_gene_name", "product", "number_of_genomes_present", "prevalence_percentage", "frequency_class", "graph_degree"]], ctx.output_dir / "02_tables/top_named_gene_clusters.tsv")
    product_col = master["product"].fillna("").replace("", np.nan)
    top_products = product_col.value_counts().head(200).reset_index()
    top_products.columns = ["product", "gene_cluster_count"]
    write_tsv(top_products, ctx.output_dir / "02_tables/top_products.tsv")
    conserved_hyp = master[(master["frequency_class"].isin(["core", "soft_core"])) & (master["hypothetical_or_uncharacterized"].fillna(False))].sort_values(["prevalence_percentage", "graph_degree"], ascending=False)
    write_tsv(conserved_hyp, ctx.output_dir / "02_tables/conserved_hypothetical_clusters.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/conserved_hypothetical_clusters.tsv", "Table_S8_conserved_hypothetical_clusters.tsv")
    high_deg_un = master[master["hypothetical_or_uncharacterized"].fillna(False)].sort_values("graph_degree", ascending=False).head(500)
    write_tsv(high_deg_un, ctx.output_dir / "02_tables/high_degree_uncharacterized_clusters.tsv")
    length_cols = [c for c in ["sequence_length", "protein_length", "dna_length_min", "dna_length_max", "dna_length_mean", "graph_length_mean"] if c in master.columns]
    length_summary = master.groupby("frequency_class", observed=True)[length_cols].agg(["count", "median", "mean", "min", "max"]).reset_index()
    length_summary.columns = ["_".join(map(str, c)).strip("_") for c in length_summary.columns]
    write_tsv(length_summary, ctx.output_dir / "02_tables/gene_length_summary.tsv")
    annot_summary = master.groupby("frequency_class", observed=True).agg(
        gene_clusters=("gene_cluster_id", "count"),
        uncharacterized=("hypothetical_or_uncharacterized", "sum"),
        median_prevalence=("number_of_genomes_present", "median"),
        median_graph_degree=("graph_degree", "median"),
    ).reset_index()
    annot_summary["uncharacterized_percentage"] = annot_summary["uncharacterized"] / annot_summary["gene_clusters"] * 100
    write_tsv(annot_summary, ctx.output_dir / "02_tables/annotation_summary.tsv")
    setup_plot()
    plot_annotation_figures(ctx, master, annot_summary, top_products, conserved_hyp)


def plot_annotation_figures(ctx: Context, master: pd.DataFrame, annot_summary: pd.DataFrame, top_products: pd.DataFrame, conserved_hyp: pd.DataFrame) -> None:
    source_data(ctx, "annotation_landscape_source_data.tsv", annot_summary)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    vals = annot_summary.set_index("frequency_class").reindex(FREQ_ORDER)
    ax.bar(vals.index, 100 - vals["uncharacterized_percentage"], label="Named/characterized", color="#0072B2")
    ax.bar(vals.index, vals["uncharacterized_percentage"], bottom=100 - vals["uncharacterized_percentage"], label="Hypothetical/uncharacterized", color="#7F7F7F")
    ax.set_ylabel("Gene clusters (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Annotation status")
    ax.legend(frameon=False)
    save_figure(ctx, fig, "annotation_status_by_frequency_class")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    data = master["sequence_length"].dropna()
    ax.hist(data, bins=60, color="#009E73")
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Gene clusters")
    ax.set_title("Gene length distribution")
    save_figure(ctx, fig, "gene_protein_length_distribution")

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    sample = master.dropna(subset=["sequence_length"]).sample(min(8000, master["sequence_length"].notna().sum()), random_state=42)
    ax.scatter(sample["sequence_length"], sample["prevalence_percentage"], s=4, alpha=0.25, color="#0072B2")
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Prevalence (%)")
    ax.set_title("Prevalence versus length")
    save_figure(ctx, fig, "prevalence_vs_sequence_length")

    top = top_products.head(15).iloc[::-1]
    source_data(ctx, "top_products_source_data.tsv", top_products)
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.barh(top["product"].map(lambda x: safe_label(x, 45)), top["gene_cluster_count"], color="#E69F00")
    ax.set_xlabel("Gene clusters")
    ax.set_title("Top recurrent product annotations")
    save_figure(ctx, fig, "top_recurrent_named_functions")

    top_hyp = conserved_hyp.head(25)
    source_data(ctx, "conserved_uncharacterized_source_data.tsv", top_hyp)
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.scatter(top_hyp["prevalence_percentage"], top_hyp["graph_degree"].fillna(0), s=18, color="#CC79A7")
    ax.set_xlabel("Prevalence (%)")
    ax.set_ylabel("Graph degree")
    ax.set_title("Conserved uncharacterized clusters")
    save_figure(ctx, fig, "conserved_uncharacterized_overview")


def accessory_genome_structure(ctx: Context) -> None:
    master = load_master(ctx)
    genes, genomes, matrix = load_rtab(ctx)
    prev = matrix.sum(axis=1)
    n = len(genomes)
    filt = (prev > 1) & (prev < 0.99 * n)
    filtered_genes = np.asarray(genes)[filt]
    X = sparse.csr_matrix(matrix[filt, :].T.astype(np.uint8))
    log(ctx, "accessory", f"Filtered accessory matrix: {X.shape[0]} genomes x {X.shape[1]} genes")
    intersections = (X @ X.T).toarray().astype(float)
    sums = np.asarray(X.sum(axis=1)).ravel().astype(float)
    unions = sums[:, None] + sums[None, :] - intersections
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = 1.0 - intersections / unions
    jaccard[~np.isfinite(jaccard)] = 0
    np.fill_diagonal(jaccard, 0)
    coords = pcoa(jaccard, n_components=5)
    coord_df = pd.DataFrame(coords, columns=[f"PCoA{i}" for i in range(1, coords.shape[1] + 1)])
    coord_df.insert(0, "genome_id", genomes)
    try:
        u, s, vt = svds(X.astype(float), k=5)
        order = np.argsort(s)[::-1]
        svd_coords = u[:, order] * s[order]
        for i in range(svd_coords.shape[1]):
            coord_df[f"SVD{i+1}"] = svd_coords[:, i]
    except Exception as exc:
        for i in range(1, 6):
            coord_df[f"SVD{i}"] = np.nan
        write_skip(ctx, "02_tables/accessory_svd.skip", f"Truncated SVD failed: {exc}")
    write_tsv(coord_df, ctx.output_dir / "02_tables/accessory_ordination_coordinates.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/accessory_ordination_coordinates.tsv", "Table_S14_accessory_genome_ordination_coordinates.tsv")

    condensed = squareform(jaccard, checks=False)
    Z = hierarchy.linkage(condensed, method="average")
    stabilities = []
    cluster_df = pd.DataFrame({"genome_id": genomes})
    for k in range(2, 9):
        labels = hierarchy.fcluster(Z, t=k, criterion="maxclust")
        sil = silhouette_from_distance(jaccard, labels)
        stabilities.append({"k": k, "mean_silhouette": sil, "method": "average_linkage_jaccard"})
        cluster_df[f"cluster_k{k}"] = labels
    best_k = max(stabilities, key=lambda r: r["mean_silhouette"])["k"]
    cluster_df["selected_cluster"] = cluster_df[f"cluster_k{best_k}"]
    write_tsv(cluster_df, ctx.output_dir / "02_tables/accessory_genome_clusters.tsv")
    write_tsv(pd.DataFrame(stabilities), ctx.output_dir / "02_tables/accessory_cluster_stability.tsv")
    dist_summary = pd.DataFrame([{"metric": "sparse_jaccard_accessory", "genes_after_filter": int(filt.sum()), "min_distance": float(np.min(condensed)), "median_distance": float(np.median(condensed)), "mean_distance": float(np.mean(condensed)), "max_distance": float(np.max(condensed))}])
    write_tsv(dist_summary, ctx.output_dir / "02_tables/accessory_distance_summary.tsv")
    setup_plot()
    plot_accessory_figures(ctx, coord_df, cluster_df, stabilities, jaccard, Z, matrix, filtered_genes, filt, master)


def pcoa(distance: np.ndarray, n_components: int = 5) -> np.ndarray:
    n = distance.shape[0]
    D2 = distance ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 0
    vals = eigvals[positive][:n_components]
    vecs = eigvecs[:, positive][:, :n_components]
    coords = vecs * np.sqrt(vals)
    if coords.shape[1] < n_components:
        coords = np.pad(coords, ((0, 0), (0, n_components - coords.shape[1])), constant_values=0)
    return coords


def silhouette_from_distance(distance: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels)
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan
    vals = []
    for i in range(distance.shape[0]):
        same = labels == labels[i]
        other_labels = [u for u in unique if u != labels[i]]
        a = distance[i, same].sum() / max(1, same.sum() - 1)
        b = min(distance[i, labels == u].mean() for u in other_labels)
        vals.append((b - a) / max(a, b) if max(a, b) else 0)
    return float(np.mean(vals))


def plot_accessory_figures(ctx: Context, coord_df: pd.DataFrame, cluster_df: pd.DataFrame, stabilities: list[dict[str, Any]], jaccard: np.ndarray, Z: np.ndarray, matrix: np.ndarray, filtered_genes: np.ndarray, filt: np.ndarray, master: pd.DataFrame) -> None:
    source_data(ctx, "accessory_ordination_source_data.tsv", coord_df)
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.scatter(coord_df["PCoA1"], coord_df["PCoA2"], s=8, alpha=0.6, color="#0072B2")
    ax.set_xlabel("PCoA1")
    ax.set_ylabel("PCoA2")
    ax.set_title("Accessory-genome PCoA")
    save_figure(ctx, fig, "accessory_genome_pcoa")

    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.scatter(coord_df["SVD1"], coord_df["SVD2"], s=8, alpha=0.6, color="#009E73")
    ax.set_xlabel("SVD1")
    ax.set_ylabel("SVD2")
    ax.set_title("Truncated-SVD ordination")
    save_figure(ctx, fig, "accessory_genome_svd")

    prev = matrix[filt, :].mean(axis=1)
    var = prev * (1 - prev)
    top_idx = np.argsort(var)[-200:]
    heat = matrix[filt, :][top_idx, :]
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.imshow(heat, aspect="auto", interpolation="nearest", cmap="Greys", rasterized=True)
    ax.set_xlabel("Genomes")
    ax.set_ylabel("Informative accessory genes")
    ax.set_title("Accessory-gene prevalence heatmap")
    save_figure(ctx, fig, "accessory_prevalence_heatmap")
    source_data(ctx, "accessory_heatmap_gene_list.tsv", pd.DataFrame({"gene_cluster_id": filtered_genes[top_idx]}))

    per_genome = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv") if (ctx.output_dir / "02_tables/per_genome_gene_content.tsv").exists() else pd.DataFrame()
    if not per_genome.empty:
        merged = coord_df.merge(per_genome, on="genome_id", how="left")
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        sc = ax.scatter(merged["PCoA1"], merged["accessory_gene_count"], c=merged["PCoA2"], s=8, cmap="viridis", alpha=0.7)
        ax.set_xlabel("PCoA1")
        ax.set_ylabel("Accessory clusters")
        ax.set_title("Burden versus ordination")
        fig.colorbar(sc, ax=ax, label="PCoA2")
        save_figure(ctx, fig, "accessory_burden_vs_ordination")

    fig, ax = plt.subplots(figsize=(5.0, 2.8))
    hierarchy.dendrogram(Z, no_labels=True, color_threshold=None, ax=ax)
    ax.set_ylabel("Jaccard distance")
    ax.set_title("Genome clustering summary")
    save_figure(ctx, fig, "genome_clustering_dendrogram")

    stab = pd.DataFrame(stabilities)
    source_data(ctx, "accessory_cluster_stability_source_data.tsv", stab)
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.plot(stab["k"], stab["mean_silhouette"], marker="o", color="#D55E00")
    ax.set_xlabel("Clusters (k)")
    ax.set_ylabel("Mean silhouette")
    ax.set_title("Accessory-cluster stability")
    save_figure(ctx, fig, "accessory_gene_diversity_summary")


def load_or_detect_communities(ctx: Context, graph: nx.Graph) -> dict[int, int]:
    path = ctx.output_dir / "04_networks/tables/network_communities.tsv"
    if path.exists():
        df = read_tsv(path)
        if {"node_id", "community_id"}.issubset(df.columns):
            return dict(zip(df["node_id"].astype(int), df["community_id"].astype(int)))
    try:
        communities = nx.community.louvain_communities(graph, weight="weight", seed=ctx.seed)
    except Exception:
        communities = nx.community.greedy_modularity_communities(graph, weight="weight")
    mapping = {}
    for i, comm in enumerate(sorted(communities, key=len, reverse=True), start=1):
        for n in comm:
            mapping[int(n)] = i
    return mapping


def network_analysis(ctx: Context) -> None:
    node_df, edge_df, _ = parse_graph_cached(ctx, "final_graph.gml")
    pre_node_df, pre_edge_df, _ = parse_graph_cached(ctx, "pre_filt_graph.gml")
    graph = graph_from_tables(node_df, edge_df)
    n = graph.number_of_nodes()
    m = graph.number_of_edges()
    comps = sorted(nx.connected_components(graph), key=len, reverse=True)
    giant = graph.subgraph(comps[0]).copy() if comps else graph.copy()
    deg = dict(graph.degree())
    wdeg = dict(graph.degree(weight="weight"))
    core_numbers = nx.core_number(graph)
    community_map = load_or_detect_communities(ctx, graph)
    communities = defaultdict(list)
    for node, comm in community_map.items():
        communities[int(comm)].append(node)
    try:
        modularity = nx.community.modularity(graph, [set(v) for v in communities.values()], weight="weight")
    except Exception:
        modularity = np.nan
    approx_bet = nx.betweenness_centrality(graph, k=min(1000, n), seed=ctx.seed, weight=None)
    arts = set(nx.articulation_points(graph))
    bridges = list(nx.bridges(graph))
    bridge_nodes = {u for edge in bridges for u in edge}
    comp_id = {}
    comp_sizes = {}
    for i, comp in enumerate(comps, start=1):
        comp_sizes[i] = len(comp)
        for node in comp:
            comp_id[node] = i
    try:
        approx_diameter = nx.approximation.diameter(giant, seed=ctx.seed)
    except Exception:
        approx_diameter = np.nan
    try:
        approx_aspl = approximate_aspl(giant, seed=ctx.seed)
    except Exception:
        approx_aspl = np.nan
    global_stats = pd.DataFrame(
        [
            {"metric": "node_count", "value": n, "method": "exact"},
            {"metric": "edge_count", "value": m, "method": "exact"},
            {"metric": "density", "value": nx.density(graph), "method": "exact"},
            {"metric": "connected_components", "value": len(comps), "method": "exact"},
            {"metric": "giant_component_nodes", "value": len(comps[0]) if comps else 0, "method": "exact"},
            {"metric": "giant_component_proportion", "value": (len(comps[0]) / n if n else np.nan), "method": "exact"},
            {"metric": "isolated_nodes", "value": nx.number_of_isolates(graph), "method": "exact"},
            {"metric": "mean_degree", "value": statistics.mean(deg.values()) if deg else 0, "method": "exact"},
            {"metric": "median_degree", "value": statistics.median(deg.values()) if deg else 0, "method": "exact"},
            {"metric": "max_degree", "value": max(deg.values()) if deg else 0, "method": "exact"},
            {"metric": "degree_assortativity", "value": nx.degree_assortativity_coefficient(graph), "method": "exact_networkx"},
            {"metric": "average_clustering", "value": nx.average_clustering(graph), "method": "exact_networkx"},
            {"metric": "transitivity", "value": nx.transitivity(graph), "method": "exact_networkx"},
            {"metric": "maximum_k_core", "value": max(core_numbers.values()) if core_numbers else 0, "method": "exact"},
            {"metric": "approximate_giant_component_diameter", "value": approx_diameter, "method": "networkx_approximation"},
            {"metric": "approximate_average_shortest_path_length", "value": approx_aspl, "method": "sampled_shortest_paths"},
            {"metric": "articulation_points", "value": len(arts), "method": "exact"},
            {"metric": "bridges", "value": len(bridges), "method": "exact"},
            {"metric": "communities", "value": len(communities), "method": "networkx_louvain_or_greedy_fallback"},
            {"metric": "community_modularity", "value": modularity, "method": "exact_given_detected_partition"},
        ]
    )
    write_tsv(global_stats, ctx.output_dir / "04_networks/tables/network_global_statistics.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "04_networks/tables/network_global_statistics.tsv", "Table_S9_network_global_metrics.tsv")

    component_stats = pd.DataFrame([{"component_id": i, "node_count": len(comp), "edge_count": graph.subgraph(comp).number_of_edges(), "density": nx.density(graph.subgraph(comp)) if len(comp) > 1 else 0} for i, comp in enumerate(comps, start=1)])
    write_tsv(component_stats, ctx.output_dir / "04_networks/tables/network_component_statistics.tsv")
    node_metrics = pd.DataFrame(
        {
            "node_id": list(graph.nodes),
            "gene_cluster_id": [graph.nodes[n].get("gene_cluster_id", "") for n in graph.nodes],
            "degree": [deg[n] for n in graph.nodes],
            "weighted_degree": [wdeg[n] for n in graph.nodes],
            "component_id": [comp_id.get(n, pd.NA) for n in graph.nodes],
            "component_size": [comp_sizes.get(comp_id.get(n, 0), 0) for n in graph.nodes],
            "k_core_number": [core_numbers.get(n, 0) for n in graph.nodes],
            "approx_betweenness": [approx_bet.get(n, 0) for n in graph.nodes],
            "articulation_point": [n in arts for n in graph.nodes],
            "bridge_related": [n in bridge_nodes for n in graph.nodes],
            "community_id": [community_map.get(n, pd.NA) for n in graph.nodes],
        }
    )
    write_tsv(node_metrics, ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz")
    communities_df = node_metrics[["node_id", "gene_cluster_id", "community_id"]].copy()
    write_tsv(communities_df, ctx.output_dir / "04_networks/tables/network_communities.tsv")
    hubs = node_metrics.sort_values(["degree", "approx_betweenness", "k_core_number"], ascending=False).head(500)
    write_tsv(hubs, ctx.output_dir / "04_networks/tables/network_hubs.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "04_networks/tables/network_hubs.tsv", "Table_S10_network_hubs.tsv")
    write_tsv(node_metrics[node_metrics["articulation_point"]], ctx.output_dir / "04_networks/tables/network_articulation_points.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "04_networks/tables/network_articulation_points.tsv", "Table_S11_articulation_and_bridge_associated_genes.tsv")
    write_tsv(pd.DataFrame(bridges, columns=["source", "target"]), ctx.output_dir / "04_networks/tables/network_bridge_edges.tsv")
    kcore_summary = node_metrics.groupby("k_core_number").size().reset_index(name="node_count")
    write_tsv(kcore_summary, ctx.output_dir / "04_networks/tables/network_kcore_summary.tsv")
    comp = pd.DataFrame([{"graph": "final_graph", "nodes": len(node_df), "edges": len(edge_df)}, {"graph": "pre_filt_graph", "nodes": len(pre_node_df), "edges": len(pre_edge_df)}])
    write_tsv(comp, ctx.output_dir / "04_networks/tables/final_vs_prefilter_graph_comparison.tsv")
    master = load_master(ctx)
    community_characteristics(ctx, graph, node_metrics, master)
    export_networks(ctx, graph, node_metrics)


def approximate_aspl(graph: nx.Graph, seed: int = 42, samples: int = 128) -> float:
    rng = random.Random(seed)
    nodes = list(graph.nodes)
    if not nodes:
        return np.nan
    selected = rng.sample(nodes, min(samples, len(nodes)))
    values = []
    for node in selected:
        lengths = nx.single_source_shortest_path_length(graph, node)
        values.extend(lengths.values())
    return float(np.mean(values)) if values else np.nan


def community_characteristics(ctx: Context, graph: nx.Graph, node_metrics: pd.DataFrame, master: pd.DataFrame) -> None:
    nm = node_metrics.merge(master[["gene_cluster_id", "frequency_class", "hypothetical_or_uncharacterized", "product", "number_of_genomes_present", "graph_degree"]], on="gene_cluster_id", how="left")
    rows = []
    rep_rows = []
    enrich_rows = []
    for comm, sub in nm.groupby("community_id", dropna=True):
        nodes = set(sub["node_id"].astype(int))
        sg = graph.subgraph(nodes)
        freq_counts = sub["frequency_class"].value_counts().to_dict()
        central = sub.sort_values(["degree", "approx_betweenness"], ascending=False).iloc[0]
        rows.append(
            {
                "community_id": int(comm),
                "node_count": len(sub),
                "edge_count": sg.number_of_edges(),
                "internal_density": nx.density(sg) if len(sg) > 1 else 0,
                "median_prevalence": sub["number_of_genomes_present"].median(),
                "core_count": freq_counts.get("core", 0),
                "soft_core_count": freq_counts.get("soft_core", 0),
                "shell_count": freq_counts.get("shell", 0),
                "cloud_count": freq_counts.get("cloud", 0),
                "proportion_core": freq_counts.get("core", 0) / len(sub),
                "proportion_accessory": (freq_counts.get("shell", 0) + freq_counts.get("cloud", 0)) / len(sub),
                "proportion_cloud": freq_counts.get("cloud", 0) / len(sub),
                "annotation_coverage": 1 - sub["hypothetical_or_uncharacterized"].fillna(True).mean(),
                "proportion_hypothetical": sub["hypothetical_or_uncharacterized"].fillna(True).mean(),
                "central_hub_gene_cluster_id": central["gene_cluster_id"],
                "mean_k_core_number": sub["k_core_number"].mean(),
                "component_membership_mode": sub["component_id"].mode().iloc[0] if not sub["component_id"].mode().empty else pd.NA,
            }
        )
        reps = sub[~sub["hypothetical_or_uncharacterized"].fillna(True)].sort_values(["degree", "number_of_genomes_present"], ascending=False).head(10)
        for _, r in reps.iterrows():
            rep_rows.append({"community_id": int(comm), "gene_cluster_id": r["gene_cluster_id"], "product": r.get("product", ""), "degree": r["degree"], "prevalence": r["number_of_genomes_present"]})
        for cls in FREQ_ORDER:
            a = int((sub["frequency_class"] == cls).sum())
            b = len(sub) - a
            c = int((nm["frequency_class"] == cls).sum()) - a
            d = len(nm) - len(sub) - c
            if min(a + b, c + d, a + c, b + d) > 0:
                try:
                    odds, p = stats.fisher_exact([[a, b], [c, d]])
                except Exception:
                    odds, p = np.nan, np.nan
                enrich_rows.append({"community_id": int(comm), "category": cls, "test": "Fisher exact one category versus remaining nodes", "a_in_community": a, "b_other_in_community": b, "c_in_background": c, "d_other_background": d, "odds_ratio": odds, "raw_p_value": p})
    char = pd.DataFrame(rows).sort_values("node_count", ascending=False)
    reps = pd.DataFrame(rep_rows)
    enrich = pd.DataFrame(enrich_rows)
    if not enrich.empty:
        enrich["fdr_adjusted_p_value"] = benjamini_hochberg(enrich["raw_p_value"])
        enrich["comparisons"] = len(enrich)
    write_tsv(char, ctx.output_dir / "04_networks/tables/community_characteristics.tsv")
    write_tsv(enrich, ctx.output_dir / "04_networks/tables/community_enrichment_results.tsv")
    write_tsv(reps, ctx.output_dir / "04_networks/tables/community_representative_genes.tsv")
    summary = char.merge(reps.groupby("community_id")["product"].apply(lambda s: "; ".join(map(str, s.head(5)))).reset_index(name="representative_named_products"), on="community_id", how="left")
    write_tsv(summary, ctx.output_dir / "04_networks/tables/community_annotation_summary.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "04_networks/tables/community_characteristics.tsv", "Table_S12_network_communities.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "04_networks/tables/community_enrichment_results.tsv", "Table_S13_community_enrichment_results.tsv")


def export_networks(ctx: Context, graph: nx.Graph, node_metrics: pd.DataFrame) -> None:
    # Attach scalar metrics for export.
    metrics = node_metrics.set_index("node_id").to_dict("index")
    for n in graph.nodes:
        for k, v in metrics.get(n, {}).items():
            if pd.isna(v) if not isinstance(v, (list, dict, tuple)) else False:
                continue
            graph.nodes[n][k] = str(v) if isinstance(v, pd._libs.missing.NAType) else v
    nx.write_graphml(graph, ctx.output_dir / "04_networks/graph_exports/full_annotated_network.graphml")
    write_tsv(node_metrics, ctx.output_dir / "04_networks/cytoscape/cytoscape_node_table.tsv")
    edge_table = pd.DataFrame([(u, v, d.get("weight", 1)) for u, v, d in graph.edges(data=True)], columns=["source", "target", "weight"])
    write_tsv(edge_table, ctx.output_dir / "04_networks/cytoscape/cytoscape_edge_table.tsv")
    save_text(
        """
        # Cytoscape Style Guide

        Suggested mappings:

        - Node label: `gene_cluster_id`
        - Node size: log-scaled `degree` or `weighted_degree`
        - Node color: `frequency_class` when merged with `gene_cluster_master.tsv.gz`, or `community_id`
        - Edge width: `weight`
        - Filter full graph by `k_core_number`, `degree`, or selected communities before manual layout.
        """,
        ctx.output_dir / "04_networks/cytoscape/cytoscape_style_guide.md",
    )
    max_k = int(node_metrics["k_core_number"].max())
    high_nodes = set(node_metrics[node_metrics["k_core_number"] >= max(1, max_k - 1)]["node_id"].astype(int))
    high = graph.subgraph(high_nodes).copy()
    nx.write_graphml(high, ctx.output_dir / "04_networks/graph_exports/high_kcore_network.graphml")
    nx.write_gml(high, ctx.output_dir / "04_networks/graph_exports/high_kcore_network.gml")
    top_nodes = set(node_metrics.sort_values(["degree", "approx_betweenness"], ascending=False).head(300)["node_id"].astype(int))
    bridge_nodes = set(node_metrics[node_metrics["bridge_related"]]["node_id"].astype(int))
    hb_nodes = top_nodes | set(list(bridge_nodes)[:300])
    hb = graph.subgraph(hb_nodes).copy()
    nx.write_graphml(hb, ctx.output_dir / "04_networks/graph_exports/hub_bridge_network.graphml")
    nx.write_gml(hb, ctx.output_dir / "04_networks/graph_exports/hub_bridge_network.gml")
    nx.node_link_data(hb, edges="edges")
    save_json(nx.node_link_data(hb, edges="edges"), ctx.output_dir / "04_networks/graph_exports/hub_bridge_network.json")
    backbone_edges = []
    comm_centers = node_metrics.sort_values("degree", ascending=False).groupby("community_id").head(1)
    center_nodes = set(comm_centers["node_id"].astype(int))
    for u, v, d in graph.edges(data=True):
        if u in center_nodes and v in center_nodes:
            backbone_edges.append((u, v, d))
    backbone = nx.Graph()
    backbone.add_nodes_from((n, graph.nodes[n]) for n in center_nodes)
    backbone.add_edges_from(backbone_edges)
    nx.write_graphml(backbone, ctx.output_dir / "04_networks/graph_exports/community_backbone.graphml")


def network_visualization(ctx: Context) -> None:
    master = load_master(ctx)
    node_metrics = read_tsv(ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz")
    component_stats = read_tsv(ctx.output_dir / "04_networks/tables/network_component_statistics.tsv")
    global_stats = read_tsv(ctx.output_dir / "04_networks/tables/network_global_statistics.tsv")
    kcore = read_tsv(ctx.output_dir / "04_networks/tables/network_kcore_summary.tsv")
    community = read_tsv(ctx.output_dir / "04_networks/tables/community_characteristics.tsv")
    setup_plot()
    source_data(ctx, "network_metrics_source_data.tsv", node_metrics[["node_id", "gene_cluster_id", "degree", "weighted_degree", "k_core_number", "approx_betweenness", "community_id"]])
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.hist2d(node_metrics["degree"], node_metrics["weighted_degree"], bins=60, cmap="viridis")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Weighted degree")
    ax.set_title("Full-network density representation")
    save_figure(ctx, fig, "full_network_density_representation")
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(np.arange(1, len(component_stats) + 1), component_stats["node_count"], color="#0072B2")
    ax.set_yscale("log")
    ax.set_xlabel("Component rank")
    ax.set_ylabel("Nodes (log scale)")
    ax.set_title("Connected-component sizes")
    save_figure(ctx, fig, "connected_component_size_distribution")
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    deg_counts = node_metrics["degree"].value_counts().sort_index()
    ax.bar(deg_counts.index, deg_counts.values, color="#009E73")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Nodes")
    ax.set_title("Degree distribution")
    save_figure(ctx, fig, "degree_distribution_linear")
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.scatter(deg_counts.index, deg_counts.values, s=12, color="#009E73")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Degree (log)")
    ax.set_ylabel("Nodes (log)")
    ax.set_title("Degree distribution")
    save_figure(ctx, fig, "degree_distribution_loglog")
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.bar(kcore["k_core_number"], kcore["node_count"], color="#E69F00")
    ax.set_xlabel("k-core number")
    ax.set_ylabel("Nodes")
    ax.set_title("k-core distribution")
    save_figure(ctx, fig, "k_core_distribution")
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(np.arange(1, len(community) + 1), community.sort_values("node_count", ascending=False)["node_count"], color="#CC79A7")
    ax.set_yscale("log")
    ax.set_xlabel("Community rank")
    ax.set_ylabel("Nodes (log scale)")
    ax.set_title("Community-size distribution")
    save_figure(ctx, fig, "community_size_distribution")
    plot_selected_networks(ctx, node_metrics, master)
    plot_module_heatmap(ctx, community)
    plot_network_composites(ctx)


def plot_selected_networks(ctx: Context, node_metrics: pd.DataFrame, master: pd.DataFrame) -> None:
    node_df, edge_df, _ = parse_graph_cached(ctx, "final_graph.gml")
    graph = graph_from_tables(node_df, edge_df)
    merged = node_metrics.merge(master[["gene_cluster_id", "frequency_class", "hypothetical_or_uncharacterized"]], on="gene_cluster_id", how="left")
    selections = {
        "high_kcore_backbone": set(merged[merged["k_core_number"] >= max(merged["k_core_number"].max() - 1, 1)]["node_id"].astype(int)),
        "hub_bridge_subnetwork": set(merged.sort_values(["degree", "approx_betweenness"], ascending=False).head(250)["node_id"].astype(int)) | set(merged[merged["bridge_related"]].head(250)["node_id"].astype(int)),
        "conserved_core_neighborhood": set(merged[merged["frequency_class"].isin(["core", "soft_core"])].sort_values("degree", ascending=False).head(250)["node_id"].astype(int)),
        "accessory_gene_neighborhood": set(merged[merged["frequency_class"].isin(["shell", "cloud"])].sort_values("degree", ascending=False).head(250)["node_id"].astype(int)),
        "uncharacterized_high_centrality": set(merged[merged["hypothetical_or_uncharacterized"].fillna(False)].sort_values(["approx_betweenness", "degree"], ascending=False).head(250)["node_id"].astype(int)),
    }
    for name, nodes in selections.items():
        expanded = set(nodes)
        for n in list(nodes):
            expanded.update(list(graph.neighbors(n))[:5])
        sg = graph.subgraph(expanded).copy()
        if sg.number_of_nodes() > 600:
            keep = sorted(sg.nodes, key=lambda x: graph.degree(x), reverse=True)[:600]
            sg = sg.subgraph(keep).copy()
        draw_subgraph(ctx, sg, name)
        if Network and sg.number_of_nodes() <= 600:
            html_path = ctx.output_dir / "06_interactive" / f"{name}.html"
            net = Network(height="760px", width="100%", bgcolor="#ffffff", font_color="#222222")
            net.from_nx(sg)
            net.write_html(str(html_path), notebook=False, open_browser=False)


def draw_subgraph(ctx: Context, sg: nx.Graph, basename: str) -> None:
    if sg.number_of_nodes() == 0:
        write_skip(ctx, f"03_figures/individual/{basename}.pdf", "No nodes selected for this subnetwork.")
        return
    pos_path = ctx.output_dir / "04_networks/graph_exports" / f"{basename}_layout.tsv"
    if pos_path.exists():
        pos_df = read_tsv(pos_path)
        pos = {int(r.node_id): (r.x, r.y) for r in pos_df.itertuples()}
    else:
        pos = nx.spring_layout(sg, seed=ctx.seed, iterations=80, weight="weight")
        write_tsv(pd.DataFrame([{"node_id": n, "x": xy[0], "y": xy[1]} for n, xy in pos.items()]), pos_path)
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    degrees = dict(sg.degree())
    sizes = [10 + 5 * math.log1p(degrees[n]) for n in sg.nodes]
    nx.draw_networkx_edges(sg, pos, ax=ax, width=0.3, alpha=0.25, edge_color="#555555")
    nx.draw_networkx_nodes(sg, pos, ax=ax, node_size=sizes, node_color="#0072B2", alpha=0.85, linewidths=0)
    labels = {n: safe_label(sg.nodes[n].get("gene_cluster_id", n), 20) for n in sorted(sg.nodes, key=lambda x: degrees[x], reverse=True)[:20]}
    nx.draw_networkx_labels(sg, pos, labels=labels, font_size=5, ax=ax)
    ax.set_title(basename.replace("_", " "))
    ax.set_axis_off()
    save_figure(ctx, fig, basename)
    write_tsv(pd.DataFrame([{"source": u, "target": v, "weight": d.get("weight", 1)} for u, v, d in sg.edges(data=True)]), ctx.output_dir / "04_networks/graph_exports" / f"{basename}_edges.tsv")
    write_tsv(pd.DataFrame([{"node_id": n, **sg.nodes[n]} for n in sg.nodes]), ctx.output_dir / "04_networks/graph_exports" / f"{basename}_nodes.tsv")


def plot_module_heatmap(ctx: Context, community: pd.DataFrame) -> None:
    top = community.sort_values("node_count", ascending=False).head(30).copy()
    cols = ["proportion_core", "proportion_accessory", "proportion_cloud", "annotation_coverage", "proportion_hypothetical", "mean_k_core_number"]
    data = top[cols].to_numpy(float).copy()
    data[:, -1] = data[:, -1] / np.nanmax(data[:, -1]) if np.nanmax(data[:, -1]) else data[:, -1]
    source_data(ctx, "module_summary_heatmap_source_data.tsv", top[["community_id"] + cols])
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    im = ax.imshow(data, aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["community_id"])
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right")
    ax.set_ylabel("Community")
    ax.set_title("Module summary heatmap")
    fig.colorbar(im, ax=ax, label="Scaled value")
    save_figure(ctx, fig, "module_summary_heatmap")


def plot_network_composites(ctx: Context) -> None:
    # Final composite assembly happens in generate_publication_figures.
    return


def sparse_cooccurrence(ctx: Context, max_genes: int = 1200, min_prevalence: int = 20, max_prev_frac: float = 0.95, min_shared: int = 10) -> None:
    master = load_master(ctx)
    genes, genomes, matrix = load_rtab(ctx)
    n = len(genomes)
    prev = matrix.sum(axis=1)
    max_prev = int(max_prev_frac * n)
    mask = (prev >= min_prevalence) & (prev <= max_prev)
    variable_score = (prev / n) * (1 - prev / n)
    idx = np.where(mask)[0]
    if len(idx) > max_genes:
        idx = idx[np.argsort(variable_score[idx])[-max_genes:]]
    selected = matrix[idx, :].astype(np.uint8)
    selected_genes = np.asarray(genes)[idx]
    filter_summary = pd.DataFrame(
        [
            {"stage": "all_genes", "genes_remaining": len(genes), "pairs_possible": len(genes) * (len(genes) - 1) // 2},
            {"stage": "prevalence_filter", "genes_remaining": int(mask.sum()), "pairs_possible": int(mask.sum()) * (int(mask.sum()) - 1) // 2},
            {"stage": "variance_cap_for_pair_testing", "genes_remaining": len(idx), "pairs_possible": len(idx) * (len(idx) - 1) // 2},
        ]
    )
    write_tsv(filter_summary, ctx.output_dir / "02_tables/cooccurrence_filtering_summary.tsv")
    X = selected.astype(bool)
    intersections = X @ X.T
    prevalences = X.sum(axis=1)
    rows = []
    total_pairs = len(idx) * (len(idx) - 1) // 2
    for a in range(len(idx)):
        inter_row = intersections[a, a + 1 :]
        candidates = np.where(inter_row >= min_shared)[0] + a + 1
        # Add also low-overlap candidates among a bounded high-variance set for mutual-exclusion screening.
        low_candidates = np.where(inter_row <= max(1, min_shared // 3))[0][:50] + a + 1
        for b in np.unique(np.concatenate([candidates, low_candidates])):
            n11 = int(intersections[a, b])
            n10 = int(prevalences[a] - n11)
            n01 = int(prevalences[b] - n11)
            n00 = n - n11 - n10 - n01
            denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
            phi = ((n11 * n00) - (n10 * n01)) / denom if denom else np.nan
            jacc = n11 / (n11 + n10 + n01) if (n11 + n10 + n01) else np.nan
            odds = ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
            rows.append({"gene_cluster_id_a": selected_genes[a], "gene_cluster_id_b": selected_genes[b], "n11_shared": n11, "n10_a_only": n10, "n01_b_only": n01, "n00_neither": n00, "jaccard": jacc, "phi": phi, "odds_ratio_haldane_anscombe": odds})
    pairs = pd.DataFrame(rows)
    if pairs.empty:
        write_skip(ctx, "02_tables/significant_gene_cooccurrence_pairs.tsv.gz", "No pairs passed sparse co-occurrence filtering.")
        return
    pairs["association_type"] = np.where(pairs["phi"] >= 0, "positive_cooccurrence", "mutual_exclusion")
    test_candidates = pairs.reindex(pairs["phi"].abs().sort_values(ascending=False).index).head(50000).copy()
    pvals = []
    for r in test_candidates.itertuples():
        try:
            _, p = stats.fisher_exact([[r.n11_shared, r.n10_a_only], [r.n01_b_only, r.n00_neither]])
        except Exception:
            p = np.nan
        pvals.append(p)
    test_candidates["raw_p_value"] = pvals
    test_candidates["fdr_adjusted_p_value"] = benjamini_hochberg(pvals)
    test_candidates["comparisons"] = len(test_candidates)
    node_df, edge_df, _ = parse_graph_cached(ctx, "final_graph.gml")
    name_to_node = dict(zip(node_df["gene_cluster_id"], node_df["node_id"]))
    neighborhood_pairs = {tuple(sorted((int(u), int(v)))) for u, v in zip(edge_df["source"], edge_df["target"])}
    def has_neighborhood(a: str, b: str) -> bool:
        if a not in name_to_node or b not in name_to_node:
            return False
        return tuple(sorted((int(name_to_node[a]), int(name_to_node[b])))) in neighborhood_pairs
    test_candidates["also_panaroo_neighborhood_edge"] = [has_neighborhood(a, b) for a, b in zip(test_candidates["gene_cluster_id_a"], test_candidates["gene_cluster_id_b"])]
    sig = test_candidates[test_candidates["fdr_adjusted_p_value"] <= 0.05].sort_values(["association_type", "fdr_adjusted_p_value", "phi"])
    write_tsv(sig, ctx.output_dir / "02_tables/significant_gene_cooccurrence_pairs.tsv.gz")
    write_tsv(test_candidates[test_candidates["phi"] > 0].sort_values(["fdr_adjusted_p_value", "phi"], ascending=[True, False]).head(500), ctx.output_dir / "02_tables/top_positive_associations.tsv")
    write_tsv(test_candidates[test_candidates["phi"] < 0].sort_values(["fdr_adjusted_p_value", "phi"], ascending=[True, True]).head(500), ctx.output_dir / "02_tables/top_negative_associations.tsv")
    overlap = pd.DataFrame(
        [
            {"metric": "candidate_pairs_tested", "value": len(test_candidates)},
            {"metric": "significant_pairs_fdr_0.05", "value": len(sig)},
            {"metric": "significant_pairs_also_neighborhood_edges", "value": int(sig["also_panaroo_neighborhood_edge"].sum()) if not sig.empty else 0},
            {"metric": "total_possible_pairs_selected_genes", "value": total_pairs},
        ]
    )
    write_tsv(overlap, ctx.output_dir / "02_tables/cooccurrence_vs_neighborhood_overlap.tsv")
    copy_for_supplement(ctx, ctx.output_dir / "02_tables/significant_gene_cooccurrence_pairs.tsv.gz", "Table_S15_filtered_gene_cooccurrence_associations.tsv.gz")


def build_knowledge_graph(ctx: Context) -> None:
    master = load_master(ctx)
    node_metrics = read_tsv(ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz") if (ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz").exists() else pd.DataFrame()
    nodes = []
    edges = []
    def add_node(node_id: str, typ: str, source: str, evidence: str, **attrs: Any) -> None:
        d = {"id": node_id, "type": typ, "source_file": source, "evidence": evidence, "provenance": "publication_deep_analysis"}
        d.update(attrs)
        nodes.append(d)
    def add_edge(edge_id: str, source: str, target: str, typ: str, source_file: str, evidence: str, **attrs: Any) -> None:
        d = {"id": edge_id, "source": source, "target": target, "type": typ, "source_file": source_file, "evidence": evidence, "provenance": "publication_deep_analysis"}
        d.update(attrs)
        edges.append(d)
    for r in master.itertuples():
        gid = f"GeneCluster:{r.gene_cluster_id}"
        add_node(gid, "GeneCluster", "gene_cluster_master.tsv.gz", "one row per Rtab gene cluster", prevalence=getattr(r, "number_of_genomes_present", None), frequency_class=getattr(r, "frequency_class", None), graph_degree=getattr(r, "graph_degree", None))
        fc = f"FrequencyClass:{r.frequency_class}"
        add_edge(f"{gid}|BELONGS_TO_FREQUENCY_CLASS|{fc}", gid, fc, "BELONGS_TO_FREQUENCY_CLASS", "gene_presence_absence.Rtab", "prevalence thresholds from summary_statistics.txt")
        if pd.notna(getattr(r, "product", np.nan)) and str(r.product).strip():
            ann_id = "Annotation:" + hashlib.sha1(str(r.product).encode()).hexdigest()[:16]
            add_node(ann_id, "Annotation", "gene_presence_absence_roary.csv/final_graph.gml", "product or graph description text", label=str(r.product)[:250])
            add_edge(f"{gid}|ANNOTATED_AS|{ann_id}", gid, ann_id, "ANNOTATED_AS", "gene_presence_absence_roary.csv/final_graph.gml", "annotation/product text")
        if pd.notna(getattr(r, "community_id", np.nan)):
            cid = f"NetworkCommunity:{int(r.community_id)}"
            add_edge(f"{gid}|MEMBER_OF_COMMUNITY|{cid}", gid, cid, "MEMBER_OF_COMMUNITY", "final_graph.gml", "network community detection")
        if pd.notna(getattr(r, "connected_component", np.nan)):
            cc = f"ConnectedComponent:{int(r.connected_component)}"
            add_edge(f"{gid}|MEMBER_OF_COMPONENT|{cc}", gid, cc, "MEMBER_OF_COMPONENT", "final_graph.gml", "connected-component analysis")
    for cls in FREQ_ORDER:
        add_node(f"FrequencyClass:{cls}", "FrequencyClass", "summary_statistics.txt", "Panaroo threshold class", label=cls)
    if not node_metrics.empty:
        for comm in sorted(node_metrics["community_id"].dropna().unique()):
            add_node(f"NetworkCommunity:{int(comm)}", "NetworkCommunity", "network_communities.tsv", "detected community", label=f"Community {int(comm)}")
        for comp in sorted(node_metrics["component_id"].dropna().unique()):
            add_node(f"ConnectedComponent:{int(comp)}", "ConnectedComponent", "network_node_metrics.tsv.gz", "connected component", label=f"Component {int(comp)}")
    genes, genomes, matrix = load_rtab(ctx)
    for genome in genomes:
        add_node(f"Genome:{genome}", "Genome", "gene_presence_absence.Rtab", "Rtab genome column", label=genome)
    # Compact presence edges for representative genomes only.
    per_genome = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv") if (ctx.output_dir / "02_tables/per_genome_gene_content.tsv").exists() else pd.DataFrame({"genome_id": genomes})
    selected_genomes = list(per_genome.sort_values("accessory_gene_count", ascending=False).head(25)["genome_id"]) if "accessory_gene_count" in per_genome.columns else []
    selected_genomes += list(per_genome.sort_values("accessory_gene_count", ascending=True).head(25)["genome_id"]) if "accessory_gene_count" in per_genome.columns else []
    selected_genomes = list(dict.fromkeys(selected_genomes + genomes[:50]))
    genome_idx = {g: i for i, g in enumerate(genomes)}
    for genome in selected_genomes:
        j = genome_idx.get(genome)
        if j is None:
            continue
        present = np.asarray(genes)[matrix[:, j] == 1]
        for gene in present:
            add_edge(f"GeneCluster:{gene}|PRESENT_IN_GENOME|Genome:{genome}", f"GeneCluster:{gene}", f"Genome:{genome}", "PRESENT_IN_GENOME", "gene_presence_absence.Rtab", "selected representative genome presence edge")
    node_df, edge_df, _ = parse_graph_cached(ctx, "final_graph.gml")
    id_to_name = dict(zip(node_df["node_id"].astype(int), node_df["gene_cluster_id"].astype(str)))
    for r in edge_df.itertuples():
        a = id_to_name.get(int(r.source))
        b = id_to_name.get(int(r.target))
        if a and b:
            add_edge(f"GeneCluster:{a}|NEIGHBOR_OF|GeneCluster:{b}", f"GeneCluster:{a}", f"GeneCluster:{b}", "NEIGHBOR_OF", "final_graph.gml", "Panaroo graph edge", edge_size=getattr(r, "size", np.nan))
    co_path = ctx.output_dir / "02_tables/top_positive_associations.tsv"
    if co_path.exists():
        co = read_tsv(co_path).head(5000)
        for r in co.itertuples():
            add_edge(f"GeneCluster:{r.gene_cluster_id_a}|CO_OCCURS_WITH|GeneCluster:{r.gene_cluster_id_b}", f"GeneCluster:{r.gene_cluster_id_a}", f"GeneCluster:{r.gene_cluster_id_b}", "CO_OCCURS_WITH", "top_positive_associations.tsv", "filtered sparse co-occurrence association", phi=r.phi, fdr_adjusted_p_value=getattr(r, "fdr_adjusted_p_value", np.nan))
    nodes_df = pd.DataFrame(nodes).drop_duplicates(subset=["id"])
    edges_df = pd.DataFrame(edges).drop_duplicates(subset=["id"])
    write_tsv(nodes_df, ctx.output_dir / "05_knowledge_graph/nodes/knowledge_graph_nodes.tsv.gz")
    write_tsv(edges_df, ctx.output_dir / "05_knowledge_graph/edges/knowledge_graph_edges.tsv.gz")
    write_temp_tsv(ctx, nodes_df.groupby("type").size().reset_index(name="node_count"), "Table_S16_knowledge_graph_node_summary.tsv")
    write_temp_tsv(ctx, edges_df.groupby("type").size().reset_index(name="edge_count"), "Table_S17_knowledge_graph_edge_summary.tsv")
    kg = nx.DiGraph()
    for r in nodes_df.to_dict("records"):
        kg.add_node(r["id"], **{k: str(v) for k, v in r.items() if k != "id" and pd.notna(v)})
    for r in edges_df.to_dict("records"):
        kg.add_edge(r["source"], r["target"], key=r["id"], type=r["type"])
    nx.write_graphml(kg, ctx.output_dir / "05_knowledge_graph/graph_exports/whole_pangenome_knowledge_graph.graphml")
    save_json(nx.node_link_data(kg, edges="edges"), ctx.output_dir / "05_knowledge_graph/graph_exports/whole_pangenome_knowledge_graph.json")
    save_text(knowledge_graph_schema_text(), ctx.output_dir / "05_knowledge_graph/knowledge_graph_schema.md")
    prov = pd.DataFrame([
        {"entity": "GeneCluster", "source": "gene_presence_absence.Rtab; gene_cluster_master.tsv.gz"},
        {"entity": "NEIGHBOR_OF", "source": "final_graph.gml"},
        {"entity": "CO_OCCURS_WITH", "source": "filtered sparse cooccurrence"},
        {"entity": "PRESENT_IN_GENOME", "source": "representative genome subset from gene_presence_absence.Rtab"},
    ])
    write_tsv(prov, ctx.output_dir / "05_knowledge_graph/knowledge_graph_provenance.tsv")
    build_kg_interactive(ctx, kg)
    plot_kg_schema(ctx, nodes_df, edges_df)


def write_temp_tsv(ctx: Context, df: pd.DataFrame, name: str) -> Path:
    path = ctx.output_dir / "08_supplementary" / name
    write_tsv(df, path)
    return path


def knowledge_graph_schema_text() -> str:
    return """
    # Knowledge Graph Schema

    Node types: `GeneCluster`, `Annotation`, `FrequencyClass`, `NetworkCommunity`, `ConnectedComponent`, and `Genome`.

    Edge types: `NEIGHBOR_OF`, `ANNOTATED_AS`, `BELONGS_TO_FREQUENCY_CLASS`, `MEMBER_OF_COMMUNITY`, `MEMBER_OF_COMPONENT`, `PRESENT_IN_GENOME`, and `CO_OCCURS_WITH`.

    Every node and edge includes a stable identifier, type, source file, evidence field, and provenance label. `PRESENT_IN_GENOME` edges are exported for a documented representative genome subset to keep the primary knowledge graph portable; the full presence/absence matrix remains in `gene_presence_absence.Rtab`.
    """


def build_kg_interactive(ctx: Context, kg: nx.DiGraph) -> None:
    if not Network:
        write_skip(ctx, "05_knowledge_graph/interactive/knowledge_graph_explorer.html", "pyvis is unavailable")
        return
    # Keep a small explorer around hubs and frequency classes.
    nodes = [n for n, d in kg.nodes(data=True) if d.get("type") in {"FrequencyClass", "NetworkCommunity"}][:60]
    gene_nodes = [n for n, d in kg.nodes(data=True) if d.get("type") == "GeneCluster"][:250]
    sub = kg.subgraph(nodes + gene_nodes).copy()
    net = Network(height="800px", width="100%", bgcolor="#ffffff", font_color="#222222", directed=True)
    net.from_nx(sub)
    net.write_html(str(ctx.output_dir / "05_knowledge_graph/interactive/knowledge_graph_explorer.html"), notebook=False, open_browser=False)


def plot_kg_schema(ctx: Context, nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> None:
    source_data(ctx, "knowledge_graph_summary_source_data.tsv", pd.concat([
        nodes_df.groupby("type").size().reset_index(name="count").assign(kind="node"),
        edges_df.groupby("type").size().reset_index(name="count").assign(kind="edge"),
    ], ignore_index=True))
    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    types = ["GeneCluster", "Annotation", "FrequencyClass", "NetworkCommunity", "ConnectedComponent", "Genome"]
    pos = nx.circular_layout(types)
    g = nx.DiGraph()
    g.add_nodes_from(types)
    g.add_edges_from([("GeneCluster", "Annotation"), ("GeneCluster", "FrequencyClass"), ("GeneCluster", "NetworkCommunity"), ("GeneCluster", "ConnectedComponent"), ("GeneCluster", "Genome"), ("GeneCluster", "GeneCluster")])
    nx.draw_networkx_nodes(g, pos, ax=ax, node_color="#0072B2", node_size=1000, alpha=0.9)
    nx.draw_networkx_edges(g, pos, ax=ax, arrows=True, arrowstyle="-|>", width=1.0, alpha=0.6)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=7)
    ax.set_title("Knowledge-graph architecture")
    ax.set_axis_off()
    save_figure(ctx, fig, "knowledge_graph_schema_schematic")


def generate_publication_figures(ctx: Context) -> None:
    setup_plot()
    # Build composite figures by arranging the most informative already-derived source data.
    comp = read_tsv(ctx.output_dir / "02_tables/pangenome_composition.tsv")
    occ = read_tsv(ctx.output_dir / "02_tables/occupancy_distribution.tsv")
    per = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv")
    pan = read_tsv(ctx.output_dir / "02_tables/pangenome_accumulation_curve.tsv")
    core = read_tsv(ctx.output_dir / "02_tables/core_genome_decay_curve.tsv")
    new = read_tsv(ctx.output_dir / "02_tables/new_gene_discovery_curve.tsv")
    create_figure1(ctx, comp, occ, per, pan, core)
    create_figure2(ctx, occ, new)
    create_figure3(ctx)
    create_figure4(ctx)
    create_figure5(ctx)
    create_figure6(ctx)
    create_overview(ctx, comp, occ, per, pan, core)
    figure_quality_audit(ctx)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def create_figure1(ctx: Context, comp: pd.DataFrame, occ: pd.DataFrame, per: pd.DataFrame, pan: pd.DataFrame, core: pd.DataFrame) -> None:
    source_data(ctx, "Figure1_source_data.tsv", pd.concat([
        comp.assign(panel="B"),
        occ.assign(panel="C"),
        per.assign(panel="D"),
        pan.assign(panel="E"),
        core.assign(panel="F"),
    ], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    audit = json.loads((ctx.output_dir / "00_audit/audit_summary.json").read_text())
    axs[0, 0].axis("off")
    axs[0, 0].text(0, 0.8, f"Genomes\n{audit['genomes']:,}", fontsize=12)
    axs[0, 0].text(0, 0.45, f"Gene clusters\n{audit['gene_clusters_rtab']:,}", fontsize=12)
    axs[0, 0].text(0, 0.1, f"Graph edges\n{audit['final_graph_edges']:,}", fontsize=12)
    panel_label(axs[0, 0], "A")
    main = comp[comp["category"].isin(FREQ_ORDER)]
    axs[0, 1].bar(main["category"], main["count"], color=[FREQ_COLORS[x] for x in main["category"]])
    axs[0, 1].tick_params(axis="x", rotation=30)
    axs[0, 1].set_ylabel("Gene clusters")
    panel_label(axs[0, 1], "B")
    axs[0, 2].plot(occ["genomes_present"], occ["gene_cluster_count"], color="#0072B2")
    axs[0, 2].set_xlabel("Genomes present")
    axs[0, 2].set_ylabel("Clusters")
    panel_label(axs[0, 2], "C")
    axs[1, 0].hist(per["total_detected_gene_clusters"], bins=35, color="#0072B2", alpha=0.8, label="Total")
    axs[1, 0].hist(per["accessory_gene_count"], bins=35, color="#E69F00", alpha=0.7, label="Accessory")
    axs[1, 0].set_xlabel("Clusters per genome")
    axs[1, 0].set_ylabel("Genomes")
    axs[1, 0].legend(frameon=False)
    panel_label(axs[1, 0], "D")
    axs[1, 1].plot(pan["n_genomes"], pan["pangenome_mean"], color="#0072B2")
    axs[1, 1].fill_between(pan["n_genomes"], pan["pangenome_lower95"], pan["pangenome_upper95"], color="#0072B2", alpha=0.2)
    axs[1, 1].set_xlabel("Genomes")
    axs[1, 1].set_ylabel("Pangenome")
    panel_label(axs[1, 1], "E")
    axs[1, 2].plot(core["n_genomes"], core["core_mean"], color="#009E73")
    axs[1, 2].fill_between(core["n_genomes"], core["core_lower95"], core["core_upper95"], color="#009E73", alpha=0.2)
    axs[1, 2].set_xlabel("Genomes")
    axs[1, 2].set_ylabel("Core genome")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 1. Whole-pangenome architecture", y=0.995)
    save_figure(ctx, fig, "Figure1_whole_pangenome_architecture", composite=True)


def create_figure2(ctx: Context, occ: pd.DataFrame, new: pd.DataFrame) -> None:
    heaps = read_tsv(ctx.output_dir / "02_tables/heaps_law_model_parameters.tsv")
    rare = read_tsv(ctx.output_dir / "02_tables/rare_gene_summary.tsv")
    ann = read_tsv(ctx.output_dir / "02_tables/annotation_coverage_by_frequency_class.tsv")
    hyp = read_tsv(ctx.output_dir / "02_tables/conserved_hypothetical_clusters.tsv").head(50)
    source_data(ctx, "Figure2_source_data.tsv", pd.concat([new.assign(panel="B"), rare.assign(panel="C"), occ.assign(panel="D"), ann.assign(panel="E"), hyp.assign(panel="F")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    pan_source = read_tsv(ctx.output_dir / "03_figures/figure_source_data/heaps_fit_source_data.tsv")
    axs[0, 0].plot(pan_source["n_genomes"], pan_source["observed_mean"], color="#0072B2")
    axs[0, 0].plot(pan_source["n_genomes"], pan_source["fitted"], color="black", ls="--")
    axs[0, 0].set_xlabel("Genomes")
    axs[0, 0].set_ylabel("Pangenome")
    panel_label(axs[0, 0], "A")
    axs[0, 1].plot(new["n_genomes"], new["new_gene_mean"], color="#D55E00")
    axs[0, 1].set_xlabel("Genomes")
    axs[0, 1].set_ylabel("New clusters")
    panel_label(axs[0, 1], "B")
    axs[0, 2].bar(rare["definition"], rare["count"], color="#CC79A7")
    axs[0, 2].tick_params(axis="x", rotation=35)
    axs[0, 2].set_ylabel("Clusters")
    panel_label(axs[0, 2], "C")
    rank = np.sort(occ.loc[occ.index.repeat(occ["gene_cluster_count"]), "genomes_present"].to_numpy())[::-1]
    axs[1, 0].plot(np.arange(1, len(rank) + 1), rank, color="#0072B2", lw=1)
    axs[1, 0].set_xlabel("Rank")
    axs[1, 0].set_ylabel("Genomes present")
    panel_label(axs[1, 0], "D")
    ann = ann.set_index("frequency_class").reindex(FREQ_ORDER)
    axs[1, 1].bar(ann.index, ann["annotation_coverage_percentage"], color=[FREQ_COLORS[x] for x in ann.index])
    axs[1, 1].set_ylim(0, 100)
    axs[1, 1].set_ylabel("Annotation coverage (%)")
    panel_label(axs[1, 1], "E")
    axs[1, 2].scatter(hyp["prevalence_percentage"], hyp["graph_degree"].fillna(0), s=12, color="#7F7F7F")
    axs[1, 2].set_xlabel("Prevalence (%)")
    axs[1, 2].set_ylabel("Graph degree")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 2. Pangenome openness and rare-gene landscape", y=0.995)
    save_figure(ctx, fig, "Figure2_openness_and_rare_gene_landscape", composite=True)


def create_figure3(ctx: Context) -> None:
    node_metrics = read_tsv(ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz")
    components = read_tsv(ctx.output_dir / "04_networks/tables/network_component_statistics.tsv")
    kcore = read_tsv(ctx.output_dir / "04_networks/tables/network_kcore_summary.tsv")
    source_data(ctx, "Figure3_source_data.tsv", pd.concat([node_metrics.assign(panel="A_C_D"), components.assign(panel="B"), kcore.assign(panel="D")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    axs[0, 0].hist2d(node_metrics["degree"], node_metrics["weighted_degree"], bins=50, cmap="viridis")
    axs[0, 0].set_xlabel("Degree")
    axs[0, 0].set_ylabel("Weighted degree")
    panel_label(axs[0, 0], "A")
    axs[0, 1].plot(np.arange(1, len(components) + 1), components["node_count"], color="#0072B2")
    axs[0, 1].set_yscale("log")
    axs[0, 1].set_xlabel("Component rank")
    axs[0, 1].set_ylabel("Nodes")
    panel_label(axs[0, 1], "B")
    dc = node_metrics["degree"].value_counts().sort_index()
    axs[0, 2].scatter(dc.index, dc.values, s=10, color="#009E73")
    axs[0, 2].set_xscale("log")
    axs[0, 2].set_yscale("log")
    axs[0, 2].set_xlabel("Degree")
    axs[0, 2].set_ylabel("Nodes")
    panel_label(axs[0, 2], "C")
    axs[1, 0].bar(kcore["k_core_number"], kcore["node_count"], color="#E69F00")
    axs[1, 0].set_xlabel("k-core")
    axs[1, 0].set_ylabel("Nodes")
    panel_label(axs[1, 0], "D")
    high = node_metrics.sort_values("k_core_number", ascending=False).head(200)
    axs[1, 1].scatter(high["degree"], high["approx_betweenness"], c=high["k_core_number"], cmap="viridis", s=12)
    axs[1, 1].set_xlabel("Degree")
    axs[1, 1].set_ylabel("Approx. betweenness")
    panel_label(axs[1, 1], "E")
    hb = node_metrics.sort_values(["degree", "approx_betweenness"], ascending=False).head(200)
    axs[1, 2].scatter(hb["degree"], hb["approx_betweenness"], s=12, color="#CC79A7")
    axs[1, 2].set_xlabel("Degree")
    axs[1, 2].set_ylabel("Approx. betweenness")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 3. Whole-pangenome network architecture", y=0.995)
    save_figure(ctx, fig, "Figure3_network_architecture", composite=True)


def create_figure4(ctx: Context) -> None:
    comm = read_tsv(ctx.output_dir / "04_networks/tables/community_characteristics.tsv")
    enrich = read_tsv(ctx.output_dir / "04_networks/tables/community_enrichment_results.tsv")
    hubs = read_tsv(ctx.output_dir / "04_networks/tables/network_hubs.tsv")
    master = load_master(ctx)
    source_data(ctx, "Figure4_source_data.tsv", pd.concat([comm.assign(panel="A_B"), enrich.assign(panel="D"), hubs.assign(panel="F")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    axs[0, 0].plot(np.arange(1, len(comm) + 1), comm.sort_values("node_count", ascending=False)["node_count"], color="#CC79A7")
    axs[0, 0].set_yscale("log")
    axs[0, 0].set_xlabel("Community rank")
    axs[0, 0].set_ylabel("Nodes")
    panel_label(axs[0, 0], "A")
    top = comm.sort_values("node_count", ascending=False).head(20)
    heat = top[["proportion_core", "proportion_accessory", "proportion_cloud", "annotation_coverage", "proportion_hypothetical"]].to_numpy(float)
    axs[0, 1].imshow(heat, aspect="auto", cmap="viridis")
    axs[0, 1].set_yticks(range(len(top)))
    axs[0, 1].set_yticklabels(top["community_id"])
    axs[0, 1].set_xticks(range(5))
    axs[0, 1].set_xticklabels(["core", "access.", "cloud", "annot.", "hyp."], rotation=30)
    panel_label(axs[0, 1], "B")
    axs[0, 2].scatter(top["node_count"], top["mean_k_core_number"], s=20, color="#0072B2")
    axs[0, 2].set_xlabel("Community nodes")
    axs[0, 2].set_ylabel("Mean k-core")
    panel_label(axs[0, 2], "C")
    sig = enrich.sort_values("fdr_adjusted_p_value").head(30)
    axs[1, 0].scatter(sig["odds_ratio"], -np.log10(sig["fdr_adjusted_p_value"].replace(0, np.nan)), s=12, color="#D55E00")
    axs[1, 0].set_xlabel("Odds ratio")
    axs[1, 0].set_ylabel("-log10 FDR")
    panel_label(axs[1, 0], "D")
    merged = hubs.merge(master[["gene_cluster_id", "frequency_class"]], on="gene_cluster_id", how="left")
    for cls in FREQ_ORDER:
        sub = merged[merged["frequency_class"] == cls]
        axs[1, 1].scatter(sub["degree"], sub["k_core_number"], s=12, alpha=0.6, label=cls, color=FREQ_COLORS[cls])
    axs[1, 1].set_xlabel("Degree")
    axs[1, 1].set_ylabel("k-core")
    axs[1, 1].legend(frameon=False, ncol=2)
    panel_label(axs[1, 1], "E")
    un = master[master["hypothetical_or_uncharacterized"].fillna(False)].sort_values("graph_degree", ascending=False).head(100)
    axs[1, 2].scatter(un["graph_degree"], un["prevalence_percentage"], s=12, color="#7F7F7F")
    axs[1, 2].set_xlabel("Graph degree")
    axs[1, 2].set_ylabel("Prevalence (%)")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 4. Network modules and neighborhood organization", y=0.995)
    save_figure(ctx, fig, "Figure4_network_modules", composite=True)


def create_figure5(ctx: Context) -> None:
    coord = read_tsv(ctx.output_dir / "02_tables/accessory_ordination_coordinates.tsv")
    clusters = read_tsv(ctx.output_dir / "02_tables/accessory_genome_clusters.tsv")
    per = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv")
    dist = read_tsv(ctx.output_dir / "02_tables/accessory_distance_summary.tsv")
    source_data(ctx, "Figure5_source_data.tsv", pd.concat([coord.assign(panel="A_B"), clusters.assign(panel="C_E"), per.assign(panel="D"), dist.assign(panel="F")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    axs[0, 0].scatter(coord["PCoA1"], coord["PCoA2"], s=8, alpha=0.6, color="#0072B2")
    axs[0, 0].set_xlabel("PCoA1")
    axs[0, 0].set_ylabel("PCoA2")
    panel_label(axs[0, 0], "A")
    axs[0, 1].scatter(coord["SVD1"], coord["SVD2"], s=8, alpha=0.6, color="#009E73")
    axs[0, 1].set_xlabel("SVD1")
    axs[0, 1].set_ylabel("SVD2")
    panel_label(axs[0, 1], "B")
    axs[0, 2].hist(clusters["selected_cluster"], bins=clusters["selected_cluster"].nunique(), color="#E69F00")
    axs[0, 2].set_xlabel("Selected cluster")
    axs[0, 2].set_ylabel("Genomes")
    panel_label(axs[0, 2], "C")
    merged = coord.merge(per, on="genome_id")
    axs[1, 0].scatter(merged["PCoA1"], merged["accessory_gene_count"], s=8, alpha=0.6, color="#CC79A7")
    axs[1, 0].set_xlabel("PCoA1")
    axs[1, 0].set_ylabel("Accessory clusters")
    panel_label(axs[1, 0], "D")
    axs[1, 1].hist(per["percentage_accessory_content"], bins=40, color="#7F7F7F")
    axs[1, 1].set_xlabel("Accessory content (%)")
    axs[1, 1].set_ylabel("Genomes")
    panel_label(axs[1, 1], "E")
    axs[1, 2].bar(["median", "mean", "max"], [dist["median_distance"].iloc[0], dist["mean_distance"].iloc[0], dist["max_distance"].iloc[0]], color="#0072B2")
    axs[1, 2].set_ylabel("Jaccard distance")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 5. Accessory-genome structure", y=0.995)
    save_figure(ctx, fig, "Figure5_accessory_genome_structure", composite=True)


def create_figure6(ctx: Context) -> None:
    nodes = read_tsv(ctx.output_dir / "05_knowledge_graph/nodes/knowledge_graph_nodes.tsv.gz")
    edges = read_tsv(ctx.output_dir / "05_knowledge_graph/edges/knowledge_graph_edges.tsv.gz")
    overlap = read_tsv(ctx.output_dir / "02_tables/cooccurrence_vs_neighborhood_overlap.tsv") if (ctx.output_dir / "02_tables/cooccurrence_vs_neighborhood_overlap.tsv").exists() else pd.DataFrame({"metric": [], "value": []})
    source_data(ctx, "Figure6_source_data.tsv", pd.concat([nodes.groupby("type").size().reset_index(name="count").assign(panel="A_F"), edges.groupby("type").size().reset_index(name="count").assign(panel="B_C_D"), overlap.assign(panel="E")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8))
    node_counts = nodes.groupby("type").size().sort_values()
    axs[0, 0].barh(node_counts.index, node_counts.values, color="#0072B2")
    axs[0, 0].set_xlabel("Nodes")
    panel_label(axs[0, 0], "A")
    edge_counts = edges.groupby("type").size().sort_values()
    axs[0, 1].barh(edge_counts.index, edge_counts.values, color="#009E73")
    axs[0, 1].set_xlabel("Edges")
    panel_label(axs[0, 1], "B")
    fc = edges[edges["type"] == "BELONGS_TO_FREQUENCY_CLASS"]["target"].value_counts()
    axs[0, 2].bar(fc.index.str.replace("FrequencyClass:", ""), fc.values, color=[FREQ_COLORS.get(x.replace("FrequencyClass:", ""), "#777777") for x in fc.index])
    axs[0, 2].tick_params(axis="x", rotation=30)
    axs[0, 2].set_ylabel("Gene clusters")
    panel_label(axs[0, 2], "C")
    ann = edges[edges["type"] == "ANNOTATED_AS"].groupby("target").size().sort_values(ascending=False).head(20)
    axs[1, 0].plot(np.arange(1, len(ann) + 1), ann.values, color="#E69F00")
    axs[1, 0].set_xlabel("Annotation rank")
    axs[1, 0].set_ylabel("Gene clusters")
    panel_label(axs[1, 0], "D")
    if not overlap.empty:
        axs[1, 1].bar(overlap["metric"], overlap["value"], color="#CC79A7")
        axs[1, 1].tick_params(axis="x", rotation=35)
    panel_label(axs[1, 1], "E")
    axs[1, 2].bar(["nodes", "edges"], [len(nodes), len(edges)], color=["#0072B2", "#009E73"])
    axs[1, 2].set_ylabel("Entities")
    panel_label(axs[1, 2], "F")
    fig.suptitle("Figure 6. Pangenome knowledge graph", y=0.995)
    save_figure(ctx, fig, "Figure6_knowledge_graph", composite=True)


def create_overview(ctx: Context, comp: pd.DataFrame, occ: pd.DataFrame, per: pd.DataFrame, pan: pd.DataFrame, core: pd.DataFrame) -> None:
    source_data(ctx, "whole_pangenome_publication_overview_source_data.tsv", pd.concat([comp.assign(panel="composition"), occ.assign(panel="frequency"), per.assign(panel="per_genome"), pan.assign(panel="pangenome"), core.assign(panel="core")], ignore_index=True, sort=False))
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 4.8))
    main = comp[comp["category"].isin(FREQ_ORDER)]
    axs[0, 0].bar(main["category"], main["count"], color=[FREQ_COLORS[x] for x in main["category"]])
    axs[0, 0].tick_params(axis="x", rotation=30)
    panel_label(axs[0, 0], "A")
    axs[0, 1].plot(occ["genomes_present"], occ["gene_cluster_count"], color="#0072B2")
    panel_label(axs[0, 1], "B")
    axs[0, 2].hist(per["accessory_gene_count"], bins=35, color="#E69F00")
    panel_label(axs[0, 2], "C")
    axs[0, 3].plot(pan["n_genomes"], pan["pangenome_mean"], color="#0072B2")
    panel_label(axs[0, 3], "D")
    axs[1, 0].plot(core["n_genomes"], core["core_mean"], color="#009E73")
    panel_label(axs[1, 0], "E")
    nm = read_tsv(ctx.output_dir / "04_networks/tables/network_node_metrics.tsv.gz")
    dc = nm["degree"].value_counts().sort_index()
    axs[1, 1].scatter(dc.index, dc.values, s=8, color="#009E73")
    axs[1, 1].set_xscale("log")
    axs[1, 1].set_yscale("log")
    panel_label(axs[1, 1], "F")
    coord = read_tsv(ctx.output_dir / "02_tables/accessory_ordination_coordinates.tsv")
    axs[1, 2].scatter(coord["PCoA1"], coord["PCoA2"], s=6, alpha=0.5, color="#CC79A7")
    panel_label(axs[1, 2], "G")
    kg_nodes = read_tsv(ctx.output_dir / "05_knowledge_graph/nodes/knowledge_graph_nodes.tsv.gz")
    kg_edges = read_tsv(ctx.output_dir / "05_knowledge_graph/edges/knowledge_graph_edges.tsv.gz")
    axs[1, 3].bar(["KG nodes", "KG edges"], [len(kg_nodes), len(kg_edges)], color=["#0072B2", "#009E73"])
    panel_label(axs[1, 3], "H")
    for ax in axs.ravel():
        ax.tick_params(labelsize=6)
    fig.suptitle("Whole-pangenome publication overview", y=0.995)
    save_figure(ctx, fig, "whole_pangenome_publication_overview", composite=True)


def figure_quality_audit(ctx: Context) -> None:
    rows = []
    for path in sorted((ctx.output_dir / "03_figures").rglob("*.png")):
        if Image:
            with Image.open(path) as img:
                width, height = img.size
                dpi = img.info.get("dpi", (600, 600))
        else:
            width = height = dpi = None
        rows.append({"figure": relpath(path, ctx.output_dir), "width_px": width, "height_px": height, "dpi": dpi, "empty_plot_check": "manual_review_recommended", "clipped_text_check": "bbox_inches_tight_used", "overlap_check": "not_automated"})
    write_tsv(pd.DataFrame(rows), ctx.output_dir / "03_figures/figure_quality_audit.tsv")


def supplementary_tables(ctx: Context) -> None:
    tables = [
        ("Table_S1_input_file_inventory.tsv", "00_audit/input_file_inventory.tsv"),
        ("Table_S2_dataset_validation.tsv", "00_audit/validation_report.tsv"),
        ("Table_S3_whole_pangenome_composition.tsv", "02_tables/pangenome_composition.tsv"),
        ("Table_S4_master_gene_cluster_annotations.tsv.gz", "02_tables/gene_cluster_master.tsv.gz"),
        ("Table_S5_per_genome_gene_content_statistics.tsv", "02_tables/per_genome_gene_content.tsv"),
        ("Table_S6_gene_occupancy_distribution.tsv", "02_tables/occupancy_distribution.tsv"),
        ("Table_S7_heaps_law_model_results.tsv", "02_tables/heaps_law_model_parameters.tsv"),
        ("Table_S8_conserved_hypothetical_clusters.tsv", "02_tables/conserved_hypothetical_clusters.tsv"),
        ("Table_S9_network_global_metrics.tsv", "04_networks/tables/network_global_statistics.tsv"),
        ("Table_S10_network_hubs.tsv", "04_networks/tables/network_hubs.tsv"),
        ("Table_S11_articulation_and_bridge_associated_genes.tsv", "04_networks/tables/network_articulation_points.tsv"),
        ("Table_S12_network_communities.tsv", "04_networks/tables/community_characteristics.tsv"),
        ("Table_S13_community_enrichment_results.tsv", "04_networks/tables/community_enrichment_results.tsv"),
        ("Table_S14_accessory_genome_ordination_coordinates.tsv", "02_tables/accessory_ordination_coordinates.tsv"),
        ("Table_S15_filtered_gene_cooccurrence_associations.tsv.gz", "02_tables/significant_gene_cooccurrence_pairs.tsv.gz"),
    ]
    for out, src in tables:
        copy_for_supplement(ctx, ctx.output_dir / src, out)
    save_text(
        """
        # Supplementary Tables

        Individual TSV or TSV.GZ files are provided for Tables S1-S17 where the corresponding analysis was available. A combined Excel workbook was not produced because `openpyxl` is not installed in the available Python environment. No global package installation was attempted.
        """,
        ctx.output_dir / "08_supplementary/supplementary_table_README.md",
    )
    dd_parts = []
    for path in sorted((ctx.output_dir / "08_supplementary").glob("Table_S*.tsv*")):
        try:
            df = read_tsv(path, nrows=20)
            dd_parts.append(data_dictionary(df, path.name))
        except Exception:
            continue
    if dd_parts:
        write_tsv(pd.concat(dd_parts, ignore_index=True), ctx.output_dir / "08_supplementary/supplementary_table_data_dictionary.tsv")
    write_skip(ctx, "08_supplementary/supplementary_tables.xlsx", "openpyxl is not installed in the available Python environment.")


def manuscript_text(ctx: Context) -> None:
    audit = json.loads((ctx.output_dir / "00_audit/audit_summary.json").read_text())
    comp = read_tsv(ctx.output_dir / "02_tables/pangenome_composition.tsv").set_index("category")
    per = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv")
    heaps = read_tsv(ctx.output_dir / "02_tables/heaps_law_model_parameters.tsv")
    global_stats = read_tsv(ctx.output_dir / "04_networks/tables/network_global_statistics.tsv").set_index("metric")["value"].to_dict()
    ann = read_tsv(ctx.output_dir / "02_tables/annotation_summary.tsv")
    boot_path = ctx.output_dir / "02_tables/heaps_law_bootstrap.tsv.gz"
    boot = read_tsv(boot_path) if boot_path.exists() else pd.DataFrame()
    gamma = heaps[(heaps["model"].str.startswith("pangenome_power")) & (heaps["parameter"] == "gamma")].iloc[0]
    alpha = heaps[(heaps["model"].str.startswith("new_gene_decay")) & (heaps["parameter"] == "alpha")].iloc[0]
    def ci_text(model: str, parameter: str) -> str:
        if not boot.empty:
            vals = boot[(boot["model"] == model) & (boot["parameter"] == parameter)]["estimate"].dropna()
            if len(vals):
                return f"{vals.quantile(0.025):.4f}-{vals.quantile(0.975):.4f}"
        lo = gamma.get("bootstrap_ci95_low", np.nan) if parameter == "gamma" else alpha.get("bootstrap_ci95_low", np.nan)
        hi = gamma.get("bootstrap_ci95_high", np.nan) if parameter == "gamma" else alpha.get("bootstrap_ci95_high", np.nan)
        return f"{lo:.4f}-{hi:.4f}" if pd.notna(lo) and pd.notna(hi) else "not available"
    methods = f"""
    # Methods

    Whole-genome Panaroo pangenome outputs were analyzed under `publication_deep_analysis` without modifying the original Panaroo files. File checksums, modification times, row counts, headers, delimiters, and graph attributes were recorded before downstream analysis. Gene-cluster prevalence was recalculated from `gene_presence_absence.Rtab`.

    Frequency classes used thresholds recorded in `summary_statistics.txt`: core `(99% <= strains <= 100%)`, soft core `(95% <= strains < 99%)`, shell `(15% <= strains < 95%)`, and cloud `(0% <= strains < 15%)`. These thresholds were not inferred from external defaults.

    Randomized pangenome accumulation used {int(heaps['permutations'].dropna().iloc[0]) if 'permutations' in heaps.columns else ctx.permutations} genome-order permutations with seed {ctx.seed}. The pangenome model was `P(n) = kappa n^gamma`; the new-gene decay model was `Nnew(n) = kappa n^-alpha`. Network analyses used the post-filter `final_graph.gml`; graph adjacency is interpreted as Panaroo gene-neighborhood adjacency, not biochemical interaction.

    [CITATION REQUIRED]
    """
    results = f"""
    # Results

    ## Overall architecture of the whole P. aeruginosa pangenome

    The validated Rtab matrix contained {audit['genomes']:,} genomes and {audit['gene_clusters_rtab']:,} gene clusters. Recalculated class counts were {int(comp.loc['core','count']):,} core clusters ({comp.loc['core','percentage']:.1f}%), {int(comp.loc['soft_core','count']):,} soft-core clusters ({comp.loc['soft_core','percentage']:.1f}%), {int(comp.loc['shell','count']):,} shell clusters ({comp.loc['shell','percentage']:.1f}%), and {int(comp.loc['cloud','count']):,} cloud clusters ({comp.loc['cloud','percentage']:.1f}%) (Figure 1; Table S3).

    ## Dominance of rare and low-frequency gene clusters

    Singleton clusters numbered {int(comp.loc['singleton','count']):,} and doubletons numbered {int(comp.loc['doubleton','count']):,}. The occupancy spectrum was strongly weighted toward low-prevalence clusters, consistent with a large accessory-gene reservoir (Figure 2).

    ## Genome-to-genome variation in accessory gene content

    The median genome contained {per['total_detected_gene_clusters'].median():,.0f} detected gene clusters, including a median of {per['accessory_gene_count'].median():,.0f} accessory clusters. Outliers were reported as gene-content outliers only; no contamination or assembly-quality interpretation was made without independent QC evidence.

    ## Evidence for continuing pangenome expansion

    The fitted pangenome power-law exponent gamma was {gamma['estimate']:.4f} (bootstrap 95% CI {ci_text('pangenome_power_law', 'gamma')}). Under `P(n) = kappa n^gamma`, a positive gamma indicates continuing pangenome expansion. The new-gene decay alpha was {alpha['estimate']:.4f} (bootstrap 95% CI {ci_text('new_gene_decay', 'alpha')}). These model results support a cautious interpretation of an open pangenome for this genome collection, while recognizing that sampling density and collection composition affect rare-gene estimates.

    ## Functional distribution of conserved and accessory clusters

    Annotation coverage varied by frequency class (Table `annotation_coverage_by_frequency_class.tsv`). Clusters described as hypothetical or uncharacterized were retained as candidates for future functional investigation only when supported by annotation text; no vaccine, virulence, resistance, essentiality, or diagnostic claims were made.

    ## Gene-neighborhood network structure

    The post-filter Panaroo graph contained {int(global_stats.get('node_count', np.nan)):,} nodes and {int(global_stats.get('edge_count', np.nan)):,} edges. The giant component contained {int(global_stats.get('giant_component_nodes', np.nan)):,} nodes, and the maximum k-core was {int(global_stats.get('maximum_k_core', np.nan)):,}. Network hubs and bridge-associated genes are reported as graph-topological priorities, not direct functional interactions.

    ## Pangenome knowledge-graph construction

    A multilayer knowledge graph was built from gene clusters, annotation/product text, frequency classes, graph communities, connected components, selected genome-presence edges, Panaroo neighborhood edges, and filtered co-occurrence edges. Co-occurrence does not demonstrate causal interaction, co-regulation, horizontal transfer, or physical linkage.
    """
    discussion = """
    # Discussion and Interpretation

    The whole-genome pangenome is dominated by cloud and other low-frequency clusters, with a comparatively compact conserved core. This pattern is compatible with extensive accessory-gene diversity in the analyzed collection. The interpretation remains collection-dependent: uneven sampling, redundant lineages, assembly quality, annotation quality, and Panaroo clustering behavior can all influence the apparent abundance of rare genes.

    The Panaroo graph provides a gene-neighborhood abstraction and should not be overinterpreted as a protein-interaction or regulatory network. High-degree or high-betweenness uncharacterized nodes are useful prioritization targets for follow-up annotation or experimental work, but the present data do not establish biological mechanism.

    [CITATION REQUIRED]
    """
    limitations = """
    # Limitations

    - No clinical, environmental, geographic, phylogenetic, host, or lineage metadata were available in the Panaroo output directory; no metadata-dependent comparisons were performed.
    - Gene presence/absence can be influenced by assembly fragmentation, annotation differences, paralog handling, and pangenome clustering parameters.
    - Rare-gene estimates are sensitive to sampling density and genome collection composition.
    - Panaroo graph adjacency reflects gene-neighborhood relationships, not direct biochemical interactions.
    - Sparse co-occurrence associations are statistical associations in genome presence/absence data and do not prove functional interaction or causality.
    - Parquet and Excel workbook exports require optional packages not present in the available environment.
    """
    legends = """
    # Figure Legends

    **Figure 1. Whole-pangenome architecture.** Dataset summary, frequency-class composition, gene-frequency spectrum, per-genome gene-content distributions, pangenome accumulation, and core-genome decay.

    **Figure 2. Pangenome openness and rare-gene landscape.** Heaps' law fit, new-gene discovery, singleton/doubleton summaries, rank-prevalence structure, annotation coverage, and conserved uncharacterized clusters.

    **Figure 3. Whole-pangenome network architecture.** Network density, component sizes, degree distribution, k-core structure, high-k-core backbone metrics, and hub/bridge centrality.

    **Figure 4. Network modules and gene-neighborhood organization.** Community-size distribution, module composition, module centrality, frequency-class enrichment, conserved/accessory positions, and high-centrality uncharacterized clusters.

    **Figure 5. Accessory-genome structure.** PCoA, truncated-SVD ordination, unsupervised cluster summary, accessory burden relationships, accessory-content distribution, and Jaccard-distance summary.

    **Figure 6. Pangenome knowledge graph.** Knowledge-graph node and edge layers, frequency-class relationships, annotation links, co-occurrence/neighborhood overlap, and entity summaries.
    """
    stat = """
    # Statistical Analysis

    Pangenome accumulation used randomized genome-order permutations with empirical 95% intervals. Heaps' law parameters were estimated by non-linear least squares and bootstrap confidence intervals were obtained by resampling permutations. Community enrichment used Fisher's exact test for one category versus the remaining graph nodes with Benjamini-Hochberg FDR correction. Accessory-genome clustering used Jaccard distances and mean silhouette scores. Co-occurrence testing was restricted to a prevalence-filtered subset and used Fisher's exact tests with FDR correction after candidate filtering.
    """
    key = f"""
    # Key Findings

    - {audit['genomes']:,} genomes and {audit['gene_clusters_rtab']:,} gene clusters were independently verified from the Rtab matrix.
    - Cloud clusters dominated the pangenome ({int(comp.loc['cloud','count']):,}; {comp.loc['cloud','percentage']:.1f}%).
    - The fitted pangenome power-law gamma was {gamma['estimate']:.4f}, supporting cautious evidence for continuing pangenome expansion in this collection.
    - The final Panaroo graph contained {int(global_stats.get('node_count', np.nan)):,} nodes and {int(global_stats.get('edge_count', np.nan)):,} edges.
    """
    files = {
        "methods_whole_pangenome.md": methods,
        "results_whole_pangenome.md": results,
        "discussion_interpretation.md": discussion,
        "figure_legends.md": legends,
        "supplementary_figure_legends.md": legends.replace("# Figure Legends", "# Supplementary Figure Legends"),
        "supplementary_table_legends.md": "# Supplementary Table Legends\n\nTables S1-S17 are described in `08_supplementary/supplementary_table_README.md` and the supplementary data dictionary.\n",
        "key_findings.md": key,
        "limitations.md": limitations,
        "statistical_analysis.md": stat,
    }
    for name, text in files.items():
        save_text(text, ctx.output_dir / "07_manuscript" / name)


def html_report(ctx: Context) -> None:
    audit = json.loads((ctx.output_dir / "00_audit/audit_summary.json").read_text())
    comp = read_tsv(ctx.output_dir / "02_tables/pangenome_composition.tsv")
    final_summary_path = ctx.output_dir / "FINAL_ANALYSIS_SUMMARY.txt"
    final_summary = final_summary_path.read_text() if final_summary_path.exists() else ""
    fig_rel = [
        "03_figures/composite/Figure1_whole_pangenome_architecture.png",
        "03_figures/composite/Figure2_openness_and_rare_gene_landscape.png",
        "03_figures/composite/Figure3_network_architecture.png",
        "03_figures/composite/Figure4_network_modules.png",
        "03_figures/composite/Figure5_accessory_genome_structure.png",
        "03_figures/composite/Figure6_knowledge_graph.png",
        "03_figures/composite/whole_pangenome_publication_overview.png",
    ]
    rows_html = comp.to_html(index=False)
    template = Template(
        """
        <!doctype html>
        <html>
        <head><meta charset="utf-8"><title>Whole Pangenome Publication Report</title>
        <style>body{font-family:Arial,sans-serif;margin:2rem;line-height:1.45;color:#222} img{max-width:100%;border:1px solid #ddd;margin:1rem 0} table{border-collapse:collapse} th,td{border:1px solid #ccc;padding:.25rem .4rem} code{background:#f3f3f3;padding:.1rem .25rem}</style></head>
        <body>
        <h1>Whole P. aeruginosa Pangenome Publication Report</h1>
        <h2>Executive Summary</h2>
        <p>Validated <strong>{{ genomes }}</strong> genomes and <strong>{{ genes }}</strong> gene clusters. Final graph: <strong>{{ graph_nodes }}</strong> nodes and <strong>{{ graph_edges }}</strong> edges.</p>
        <h2>Pangenome Composition</h2>
        {{ composition }}
        <h2>Main Figures</h2>
        {% for fig in figures %}<h3>{{ fig }}</h3><img src="../{{ fig }}" alt="{{ fig }}">{% endfor %}
        <h2>Interactive Views</h2>
        <ul><li><a href="hub_bridge_subnetwork.html">Hub-bridge subnetwork</a></li><li><a href="../05_knowledge_graph/interactive/knowledge_graph_explorer.html">Knowledge graph explorer</a></li></ul>
        <h2>Output Paths</h2>
        <p>Tables: <code>02_tables/</code>; figures: <code>03_figures/</code>; networks: <code>04_networks/</code>; knowledge graph: <code>05_knowledge_graph/</code>; manuscript text: <code>07_manuscript/</code>.</p>
        <h2>Final Summary</h2><pre>{{ final_summary }}</pre>
        </body></html>
        """
    ) if Template else None
    html_text = template.render(
        genomes=f"{audit['genomes']:,}",
        genes=f"{audit['gene_clusters_rtab']:,}",
        graph_nodes=f"{audit['final_graph_nodes']:,}",
        graph_edges=f"{audit['final_graph_edges']:,}",
        composition=rows_html,
        figures=fig_rel,
        final_summary=final_summary,
    ) if template else "<html><body><h1>Report unavailable: Jinja2 missing</h1></body></html>"
    save_text(html_text, ctx.output_dir / "06_interactive/whole_pangenome_publication_report.html")


def final_summary(ctx: Context) -> None:
    audit = json.loads((ctx.output_dir / "00_audit/audit_summary.json").read_text())
    comp = read_tsv(ctx.output_dir / "02_tables/pangenome_composition.tsv").set_index("category")
    per = read_tsv(ctx.output_dir / "02_tables/per_genome_gene_content.tsv")
    heaps = read_tsv(ctx.output_dir / "02_tables/heaps_law_model_parameters.tsv")
    global_stats = read_tsv(ctx.output_dir / "04_networks/tables/network_global_statistics.tsv").set_index("metric")["value"].to_dict()
    kg_nodes = read_tsv(ctx.output_dir / "05_knowledge_graph/nodes/knowledge_graph_nodes.tsv.gz")
    kg_edges = read_tsv(ctx.output_dir / "05_knowledge_graph/edges/knowledge_graph_edges.tsv.gz")
    conserved = read_tsv(ctx.output_dir / "02_tables/conserved_hypothetical_clusters.tsv")
    validation = read_tsv(ctx.output_dir / "00_audit/validation_report.tsv")
    failed_warnings = validation[(validation["severity"] == "warning") & (validation["status"] == "fail")]
    warning_lines = "\n".join(
        f"    - {r.check}: observed {r.observed}; expected {r.expected}; {r.details if pd.notna(r.details) else ''}"
        for r in failed_warnings.itertuples()
    ) or "    - None"
    gamma = heaps[(heaps["model"].str.startswith("pangenome_power")) & (heaps["parameter"] == "gamma")].iloc[0]
    alpha = heaps[(heaps["model"].str.startswith("new_gene_decay")) & (heaps["parameter"] == "alpha")].iloc[0]
    text = f"""
    FINAL WHOLE-PANGENOME ANALYSIS SUMMARY
    =====================================
    Genomes: {audit['genomes']:,}
    Gene clusters: {audit['gene_clusters_rtab']:,}
    Core / soft-core / shell / cloud: {int(comp.loc['core','count']):,} / {int(comp.loc['soft_core','count']):,} / {int(comp.loc['shell','count']):,} / {int(comp.loc['cloud','count']):,}
    Singleton clusters: {int(comp.loc['singleton','count']):,}
    Median genes per genome: {per['total_detected_gene_clusters'].median():,.0f}
    Median accessory genes per genome: {per['accessory_gene_count'].median():,.0f}
    Heaps gamma: {gamma['estimate']:.6f}
    New-gene decay alpha: {alpha['estimate']:.6f}
    Network nodes / edges: {int(global_stats.get('node_count', np.nan)):,} / {int(global_stats.get('edge_count', np.nan)):,}
    Giant component size: {int(global_stats.get('giant_component_nodes', np.nan)):,}
    Maximum k-core: {int(global_stats.get('maximum_k_core', np.nan)):,}
    Detected communities: {int(global_stats.get('communities', np.nan)):,}
    High-priority hubs reported: {len(read_tsv(ctx.output_dir / '04_networks/tables/network_hubs.tsv')):,}
    Conserved uncharacterized clusters: {len(conserved):,}
    Knowledge-graph nodes / edges: {len(kg_nodes):,} / {len(kg_edges):,}
    Output directory: {ctx.output_dir}

    Warnings and skipped operations:
    - Validation warnings:
{warning_lines}
    - Parquet export skipped if pyarrow/fastparquet was unavailable; see .skip_reason files.
    - Combined Excel supplementary workbook skipped if openpyxl was unavailable.
    - Full gene-by-genome knowledge-graph presence edges were not exported; representative genome presence edges were exported and the complete Rtab remains the source matrix.
    - Co-occurrence associations are filtered statistical associations and do not imply causality or direct interaction.
    """
    save_text(text, ctx.output_dir / "FINAL_ANALYSIS_SUMMARY.txt")
    print(textwrap.dedent(text).strip())


def run_all(ctx: Context, args: argparse.Namespace | None = None) -> None:
    ensure_dirs(ctx)
    inventory_and_validate(ctx)
    build_master_gene_table(ctx)
    pangenome_composition(ctx)
    per_genome_gene_content(ctx)
    accumulation_and_heaps(ctx, bootstrap=getattr(args, "accumulation_bootstrap", 200) if args else 200)
    annotation_landscape(ctx)
    accessory_genome_structure(ctx)
    network_analysis(ctx)
    network_visualization(ctx)
    sparse_cooccurrence(
        ctx,
        max_genes=getattr(args, "cooccurrence_max_genes", 1200) if args else 1200,
        min_prevalence=getattr(args, "cooccurrence_min_prevalence", 20) if args else 20,
        max_prev_frac=getattr(args, "cooccurrence_max_prevalence_frac", 0.95) if args else 0.95,
        min_shared=getattr(args, "cooccurrence_min_shared", 10) if args else 10,
    )
    build_knowledge_graph(ctx)
    generate_publication_figures(ctx)
    supplementary_tables(ctx)
    manuscript_text(ctx)
    final_summary(ctx)
    html_report(ctx)
