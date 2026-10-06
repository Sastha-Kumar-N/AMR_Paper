#!/usr/bin/env python3
"""Reproducible Panaroo graph inspection, source-data export, and publication figures.

Inputs are never modified. Graph-derived measurements and annotation-derived
labels are exported as separate evidence records. Layouts do not encode distance
along a chromosome, evolutionary distance, or an experimentally verified trait.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import importlib.metadata
import itertools
import json
import logging
import math
import os
from pathlib import Path
import pickle
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Iterator

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "panaroo-figure-mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "panaroo-numba"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter
import networkx as nx
import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm

LOG = logging.getLogger("panaroo_figures")
CATEGORIES = ("core", "soft-core", "shell", "cloud")
COLORS = {"core": "#0072B2", "soft-core": "#009E73", "shell": "#E69F00", "cloud": "#CC79A7"}
EXPECTED_NODE_KEYS = {"name", "description", "annotation", "size", "members", "centroid", "seqIDs", "hasEnd", "protein", "dna", "paralog", "lengths"}
EXPECTED_EDGE_KEYS = {"size", "members", "genomeIDs"}
KEEP_NODE_KEYS = {"label", "name", "annotation", "description", "size", "paralog", "centroid", "hasEnd"}
CACHE_VERSION = 5
DEFAULTS: dict[str, Any] = {
    "total_genomes": None, "seed": 42, "organism": "Pseudomonas aeruginosa", "dataset_label": None,
    "statistics": {"enabled": True, "candidate_pairs": 50000, "fdr": 0.05,
                   "fisher_alternative": "greater", "louvain_resolution": 1.0},
    "thresholds": {"core": 0.99, "soft_core": 0.95, "shell": 0.15},
    "parser": {"nx_memory_mb": 2800, "nx_timeout_seconds": 90},
    "layout": {"algorithm": "auto", "timeout_seconds": 240},
    "figures": {"dpi": 600, "single_column_mm": 89, "double_column_mm": 183,
                "formats": ["pdf", "svg", "png", "tiff"], "histogram_bins": 50,
                "top_labels": 5, "label_genes": [], "min_edge_support": 0.05,
                "max_overview_edges": 100000, "rasterize_above_nodes": 20000,
                "focus_genes": ["mexA", "mexB", "oprM"], "neighbourhood_radius": 2,
                "max_zoom_nodes": 120, "max_zoom_labels": 12,
                "hotspot_top_n": 12, "hotspot_mad_multiplier": 3.0,
                "component_analysis": True},
}


def deep_update(base: dict, extra: dict) -> dict:
    """Merge nested configuration mappings without altering defaults."""
    result = json.loads(json.dumps(base))
    for key, value in extra.items():
        result[key] = deep_update(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else value
    return result


def read_config(path: Path | None) -> dict[str, Any]:
    """Load YAML and validate ranges before processing a graph."""
    config = deep_update(DEFAULTS, yaml.safe_load(path.read_text()) or {}) if path else deep_update(DEFAULTS, {})
    t = config["thresholds"]
    if not 0 < t["shell"] < t["soft_core"] < t["core"] <= 1:
        raise ValueError("Thresholds must satisfy 0 < shell < soft_core < core <= 1")
    f = config["figures"]
    if not 0 <= f["min_edge_support"] <= 1:
        raise ValueError("min_edge_support must be a fraction between 0 and 1")
    if f["dpi"] < 600:
        raise ValueError("Publication raster exports require dpi >= 600")
    if not {"pdf", "svg"}.issubset(f["formats"]) or not {"png", "tiff"}.intersection(f["formats"]):
        raise ValueError("formats must include pdf, svg, and at least one of png/tiff")
    if set(f["formats"]) - {"pdf", "svg", "png", "tiff"}:
        raise ValueError("Unsupported output format")
    for key in ("histogram_bins", "max_overview_edges", "max_zoom_nodes", "hotspot_top_n"):
        if f[key] < 1:
            raise ValueError(f"{key} must be positive")
    if f["neighbourhood_radius"] < 0:
        raise ValueError("neighbourhood_radius must be nonnegative")
    return config


def detect_inputs(directory: Path) -> dict[str, list[str]]:
    """Inventory current-folder inputs; do not silently search other cohorts."""
    patterns = {"gml": "*.gml", "presence_absence": "gene_presence_absence*.csv",
                "gene_data": "gene_data*.csv", "reference": "pan_genome_reference*.fa*",
                "structure": "struct_presence_absence*.Rtab"}
    return {key: [str(p.resolve()) for p in sorted(directory.glob(pattern)) if p.is_file()]
            for key, pattern in patterns.items()}


def choose_gml(inputs: dict, requested: Path | None) -> Path:
    """Select an explicit file, canonical graph, or uniquely identified clinical graph."""
    if requested:
        if not requested.is_file():
            raise FileNotFoundError(requested)
        return requested.resolve()
    files = [Path(p) for p in inputs["gml"]]
    canonical = [p for p in files if p.name == "final_graph.gml"]
    clinical = [p for p in files if "clinical" in p.stem.lower()]
    candidates = canonical or clinical or files
    if len(candidates) != 1:
        raise ValueError("Ambiguous graph selection; specify --gml explicitly")
    return candidates[0]


def sha256_file(path: Path) -> str:
    """Hash an input incrementally so caches cannot cross datasets."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tokens(path: Path) -> Iterator[str]:
    """Tokenize GML incrementally, including repeated keys and nested lists."""
    pattern = re.compile(r'"(?:[^"\\]|\\.)*"|[A-Za-z_][A-Za-z_0-9]*|[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?|[\[\]]')
    with path.open("r", encoding="ascii") as stream, tqdm(total=path.stat().st_size, unit="B", unit_scale=True, desc="Stream GML", mininterval=2) as progress:
        for line_no, line in enumerate(stream, 1):
            # Panaroo GML uses one scalar attribute per line; long sequences are
            # still a single token. Reject unexplained text rather than drop it.
            position = 0
            for match in pattern.finditer(line):
                gap = line[position:match.start()].strip()
                if gap and not gap.startswith("#"):
                    raise ValueError(f"Invalid GML at line {line_no}: {gap[:80]}")
                if gap.startswith("#"):
                    break
                yield match.group()
                position = match.end()
            tail = line[position:].strip()
            if tail and not tail.startswith("#"):
                raise ValueError(f"Unterminated or invalid GML at line {line_no}")
            progress.update(len(line))


def scalar(token: str) -> Any:
    """Decode GML scalars; HTML character references are the standard escaping."""
    if token.startswith('"'):
        return html.unescape(token[1:-1].replace('\\"', '"').replace('\\\\', '\\'))
    try:
        return int(token)
    except ValueError:
        try:
            return float(token)
        except ValueError as exc:
            raise ValueError(f"Invalid GML scalar: {token[:80]}") from exc


def read_object(stream: Iterator[str]) -> dict[str, Any]:
    """Read one bracketed record; repeated keys become lists, never last-value wins."""
    result: dict[str, Any] = {}
    repeated: set[str] = set()
    for key in stream:
        if key == "]":
            # NetworkX writes this sentinel to preserve singleton list types.
            # It is serialization syntax, never a genome ID or an annotation.
            return {k: v[1:] if isinstance(v, list) and v and v[0] == "_networkx_list_start" else v
                    for k, v in result.items()}
        if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", key):
            raise ValueError(f"Invalid GML key {key}")
        value_token = next(stream)
        value = read_object(stream) if value_token == "[" else scalar(value_token)
        if key in result:
            if key not in repeated:
                result[key] = [result[key]]
                repeated.add(key)
            result[key].append(value)
        else:
            result[key] = value
    raise ValueError("GML ended before a closing bracket")


def member_values(value: Any) -> list[str]:
    """Normalize repeated integers, semicolon strings, and list-like members."""
    if value is None:
        return []
    if isinstance(value, dict):
        return list(itertools.chain.from_iterable(member_values(v) for v in value.values()))
    if isinstance(value, (list, tuple, set)):
        return list(itertools.chain.from_iterable(member_values(v) for v in value))
    return [v.strip(" []'\"") for v in re.split(r"[;,\s]+", str(value))
            if v.strip(" []'\"") and v.strip(" []'\"") != "_networkx_list_start"]


def summarize(value: Any) -> Any:
    """Show three original attributes compactly without dumping sequences."""
    if isinstance(value, list):
        return {"type": "list", "length": len(value), "first_values": [summarize(v) for v in value[:3]]}
    if isinstance(value, dict):
        return {k: summarize(v) for k, v in value.items()}
    if isinstance(value, str) and len(value) > 160:
        return {"type": "str", "length": len(value), "prefix": value[:160]}
    return value


@dataclass
class Inspection:
    """Schema and membership provenance retained independently of compact graph."""
    node_keys: set[str] = field(default_factory=set)
    edge_keys: set[str] = field(default_factory=set)
    graph_keys: set[str] = field(default_factory=set)
    examples: list[dict] = field(default_factory=list)
    genome_ids: set[str] = field(default_factory=set)
    isolate_names: list[str] = field(default_factory=list)
    node_missing: Counter = field(default_factory=Counter)
    edge_missing: Counter = field(default_factory=Counter)
    node_size_disagreements: int = 0
    edge_size_disagreements: int = 0
    derived_node_sizes: int = 0
    derived_edge_sizes: int = 0
    parser: str = ""
    probe: str = ""


def compact_record(attributes: dict, *, is_node: bool, inspection: Inspection) -> dict:
    """Validate support and discard large sequence/membership arrays after inspection."""
    expected = EXPECTED_NODE_KEYS if is_node else EXPECTED_EDGE_KEYS
    keys = inspection.node_keys if is_node else inspection.edge_keys
    missing = inspection.node_missing if is_node else inspection.edge_missing
    keys.update(set(attributes) - {"id", "source", "target"})
    missing.update(expected - set(attributes))
    members = set(member_values(attributes.get("members", attributes.get("genomeIDs"))))
    if is_node:
        inspection.genome_ids.update(members)
        if len(inspection.examples) < 3:
            inspection.examples.append({"node_id": attributes["id"], "attributes": summarize({k: v for k, v in attributes.items() if k != "id"})})
    value = attributes.get("size")
    if value is None:
        if not members:
            raise ValueError(f"Missing size and membership support: {summarize(attributes)}")
        value = len(members)
        if is_node:
            inspection.derived_node_sizes += 1
        else:
            inspection.derived_edge_sizes += 1
    size = float(value)
    if not np.isfinite(size) or size < 0 or size != int(size):
        raise ValueError(f"Invalid size: {value}")
    if members and len(members) != size:
        if is_node:
            inspection.node_size_disagreements += 1
        else:
            inspection.edge_size_disagreements += 1
    result = {k: v for k, v in attributes.items() if k in KEEP_NODE_KEYS} if is_node else {}
    result.update(size=int(size), members_count=len(members) if members else None)
    if is_node:
        result["genome_members"] = tuple(sorted(members))
    return result


def stream_gml(path: Path) -> tuple[nx.Graph, Inspection]:
    """Parse one record at a time into a compact sparse NetworkX graph."""
    stream = tokens(path)
    if next(stream) != "graph" or next(stream) != "[":
        raise ValueError("Expected graph [ at the beginning of GML")
    graph = nx.Graph()
    inspection = Inspection(parser="streaming GML parser -> compact NetworkX Graph")
    declared: set[int] = set()
    for key in stream:
        if key == "]":
            break
        token = next(stream)
        value = read_object(stream) if token == "[" else scalar(token)
        if key == "node":
            node_id = int(value["id"])
            if node_id in declared:
                raise ValueError(f"Duplicate node id {node_id}")
            declared.add(node_id)
            graph.add_node(node_id, **compact_record(value, is_node=True, inspection=inspection))
        elif key == "edge":
            source, target = int(value["source"]), int(value["target"])
            if graph.has_edge(source, target):
                raise ValueError(f"Parallel/duplicate edge {source}-{target}; unsupported Panaroo schema")
            graph.add_edge(source, target, **compact_record(value, is_node=False, inspection=inspection))
        else:
            inspection.graph_keys.add(key)
            if key == "isolateNames":
                if value != "_networkx_list_start":
                    inspection.isolate_names.extend(member_values(value) if not isinstance(value, str) else [value])
            elif key in {"directed", "multigraph"} and int(value):
                raise ValueError("This pipeline requires an undirected, non-multigraph Panaroo graph")
    else:
        raise ValueError("Unclosed top-level graph")
    if next(stream, None) is not None:
        raise ValueError("Unexpected data after the graph")
    if set(graph) != declared:
        raise ValueError("Edges reference nodes without node records")
    return graph, inspection


def project_networkx(graph: nx.Graph) -> tuple[nx.Graph, Inspection]:
    """Preserve raw schema/examples before shrinking a successfully loaded graph."""
    if graph.is_directed() or graph.is_multigraph():
        raise ValueError("Requires an undirected simple graph")
    inspection = Inspection(parser='networkx.read_gml(path, label="id")')
    inspection.graph_keys = set(graph.graph)
    names = graph.graph.get("isolateNames", [])
    inspection.isolate_names = names if isinstance(names, list) else [names]
    for n, attrs in graph.nodes(data=True):
        compact = compact_record({"id": n, **attrs}, is_node=True, inspection=inspection)
        attrs.clear()
        attrs.update(compact)
    for u, v, attrs in graph.edges(data=True):
        compact = compact_record({"source": u, "target": v, **attrs}, is_node=False, inspection=inspection)
        attrs.clear()
        attrs.update(compact)
    graph.graph.clear()
    return graph, inspection


def networkx_probe(path: Path, destination: Path, memory_mb: int) -> None:
    """Run the requested raw reader in an isolated memory-bounded process."""
    import resource
    limit = int(memory_mb) * 1024**2
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    graph, inspection = project_networkx(nx.read_gml(path, label="id"))
    with destination.open("wb") as handle:
        pickle.dump((graph, vars(inspection)), handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_gml(path: Path, cache_dir: Path, config: dict, fingerprint: str) -> tuple[nx.Graph, Inspection]:
    """Try raw NetworkX, then streaming; reuse only locally created SHA-keyed caches."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"graph_v{CACHE_VERSION}_{fingerprint[:20]}.pickle"
    if cache.exists():
        LOG.info("Using compact graph cache %s", cache.name)
        # Only read this pipeline's local cache, never a supplied input pickle.
        with cache.open("rb") as handle:
            graph, data = pickle.load(handle)
            return graph, Inspection(**data)
    temp = cache.with_suffix(".probe")
    args = [sys.executable, str(Path(__file__).resolve()), "--_nx-probe", str(path), str(temp), str(config["parser"]["nx_memory_mb"])]
    LOG.info('Trying nx.read_gml(path, label="id") (%s MiB memory / %ss time limit)', config["parser"]["nx_memory_mb"], config["parser"]["nx_timeout_seconds"])
    try:
        probe = subprocess.run(args, capture_output=True, text=True, timeout=config["parser"]["nx_timeout_seconds"])
        if probe.returncode == 0 and temp.exists():
            with temp.open("rb") as handle:
                graph, data = pickle.load(handle)
                inspection = Inspection(**data)
            inspection.probe = "Raw NetworkX reader succeeded"
        else:
            reason = probe.stderr.strip().splitlines()[-1] if probe.stderr.strip() else f"exit {probe.returncode}"
            LOG.warning("Raw NetworkX reader failed (%s); switching to streaming parser", reason)
            graph, inspection = stream_gml(path)
            inspection.probe = f"Raw NetworkX reader failed: {reason}; streaming fallback used"
    except subprocess.TimeoutExpired:
        LOG.warning("Raw NetworkX reader timed out; switching to streaming parser")
        graph, inspection = stream_gml(path)
        inspection.probe = "Raw NetworkX reader exceeded configured timeout; streaming fallback used"
    finally:
        temp.unlink(missing_ok=True)
    with cache.open("wb") as handle:
        pickle.dump((graph, vars(inspection)), handle, protocol=pickle.HIGHEST_PROTOCOL)
    return graph, inspection


def infer_total_genomes(inspection: Inspection, configured: int | None) -> tuple[int, str]:
    """Infer cohort size from isolateNames or observed distinct member IDs, never max size."""
    header_count = len(inspection.isolate_names)
    if header_count and len(set(inspection.isolate_names)) != header_count:
        raise ValueError("Duplicate isolateNames make the genome denominator ambiguous")
    if configured is not None:
        if int(configured) != configured or configured <= 0:
            raise ValueError("total_genomes must be a positive integer")
        if header_count and int(configured) != header_count:
            raise ValueError(f"Configured total_genomes {configured} conflicts with {header_count} isolateNames")
        if len(inspection.genome_ids) > configured:
            raise ValueError("Observed genome membership exceeds configured cohort size")
        return int(configured), "config.total_genomes (validated against graph header when present)"
    if header_count:
        if len(inspection.genome_ids) > header_count:
            raise ValueError("Distinct genome members exceed isolateNames count")
        return header_count, "GML graph-level isolateNames count"
    if inspection.genome_ids:
        return len(inspection.genome_ids), "union of observed node member IDs (cannot detect isolates with no retained genes)"
    raise ValueError("Cannot infer denominator; set total_genomes in config.yaml")


def classify_prevalence(prevalence: float, thresholds: dict[str, float]) -> str:
    """Assign disjoint bins with inclusive lower boundaries at 99%, 95%, and 15%."""
    if not np.isfinite(prevalence) or not 0 <= prevalence <= 1:
        raise ValueError(f"Invalid prevalence {prevalence}")
    if prevalence >= thresholds["core"]:
        return "core"
    if prevalence >= thresholds["soft_core"]:
        return "soft-core"
    if prevalence >= thresholds["shell"]:
        return "shell"
    return "cloud"


def text_attribute(value: Any) -> str:
    """Keep all annotation alternatives; semicolon-separated names are not new genes."""
    if isinstance(value, list):
        return ";".join(map(str, value))
    return "" if value is None else str(value)


def make_node_table(graph: nx.Graph, total: int, thresholds: dict) -> pd.DataFrame:
    """Build the master node table from unfiltered topology and validated support."""
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), min(c)))
    component_id = {n: index for index, component in enumerate(components) for n in component}
    rows = []
    for n in sorted(graph):
        attrs = graph.nodes[n]
        p = attrs["size"] / total
        annotation = text_attribute(attrs.get("annotation"))
        raw_paralog = attrs.get("paralog")
        paralog = None if raw_paralog is None else str(raw_paralog).lower() in {"1", "true", "yes"}
        rows.append({"node_id": n, "panaroo_label": text_attribute(attrs.get("label")),
                     "name": text_attribute(attrs.get("name")) or f"node_{n}",
                     "annotation": annotation, "description": text_attribute(attrs.get("description")),
                     "size": attrs["size"], "prevalence": p,
                     "category": classify_prevalence(p, thresholds), "degree": graph.degree(n),
                     "paralog": paralog, "component_id": component_id[n],
                     "members_count": attrs.get("members_count"), "evidence_type": "graph_detected",
                     "label_evidence_type": "annotation_inferred"})
    return pd.DataFrame(rows).set_index("node_id", drop=False).rename_axis("_node_index")


def write_inspection(graph: nx.Graph, inspection: Inspection, total: int, denominator_source: str,
                     path: Path, fingerprint: str, destination: Path) -> None:
    """Print and save original schema, examples, membership validation, and caveats."""
    report = ["PANAROO GML INSPECTION", f"Input: {path}", f"Input SHA256: {fingerprint}",
              f"Parser: {inspection.parser}", f"Reader attempt: {inspection.probe}",
              f"Nodes: {graph.number_of_nodes():,}", f"Edges: {graph.number_of_edges():,}",
              f"Connected components (unfiltered): {nx.number_connected_components(graph):,}",
              f"Self loops: {nx.number_of_selfloops(graph):,}", f"Total genomes: {total}",
              f"Denominator source: {denominator_source}", f"Observed member IDs: {len(inspection.genome_ids)}",
              f"Graph attribute keys: {sorted(inspection.graph_keys)}",
              f"All node attribute keys: {sorted(inspection.node_keys)}",
              f"All edge attribute keys: {sorted(inspection.edge_keys)}",
              f"Expected node keys absent everywhere: {sorted(EXPECTED_NODE_KEYS - inspection.node_keys)}",
              f"Expected edge keys absent everywhere: {sorted(EXPECTED_EDGE_KEYS - inspection.edge_keys)}",
              f"Missing node attributes by count: {dict(inspection.node_missing)}",
              f"Missing edge attributes by count: {dict(inspection.edge_missing)}",
              f"Node size/member-cardinality disagreements: {inspection.node_size_disagreements}",
              f"Edge size/member-cardinality disagreements: {inspection.edge_size_disagreements}",
              f"Node/edge sizes derived from members: {inspection.derived_node_sizes}/{inspection.derived_edge_sizes}",
              "Three original example nodes (long strings/lists summarized):",
              json.dumps(inspection.examples, indent=2),
              "Memory policy: sequences are omitted from the compact cache; unique genome memberships are retained for presence/absence tests.",
              "Node size, prevalence, degree, components, and adjacency are graph-detected signals.",
              "Names, annotations, and descriptions are annotation-inferred labels; no phenotypes or insertion events are established."]
    text = "\n".join(report) + "\n"
    destination.write_text(text)
    print(text, flush=True)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--_nx-probe":
    try:
        networkx_probe(Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4]))
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)


@dataclass
class Context:
    """Shared graph, immutable master table, configuration, and provenance."""
    graph: nx.Graph
    table: pd.DataFrame
    total: int
    config: dict
    outdir: Path
    fingerprint: str
    warnings: list[str] = field(default_factory=list)
    layouts: dict[str, dict] = field(default_factory=dict)
    files: list[str] = field(default_factory=list)
    zooms: dict[str, dict] = field(default_factory=dict)
    renderers: dict[str, str] = field(default_factory=dict)
    statistics: dict[str, Any] = field(default_factory=dict)

    def warn(self, message: str) -> None:
        """Persist material limitations in logs and the machine-readable manifest."""
        if message not in self.warnings:
            self.warnings.append(message)
            LOG.warning(message)

    @property
    def options(self) -> dict:
        return self.config["figures"]


def configure_style(context: Context) -> None:
    """Use editable vector text, Okabe-Ito colors, and at least 7 pt type."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    chosen = next((f for f in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if f in installed), "DejaVu Sans")
    if chosen not in {"Arial", "Helvetica"}:
        context.warn(f"Arial/Helvetica unavailable; using metric-compatible {chosen} sans-serif.")
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [chosen],
                         "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 9,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.65, "xtick.major.width": 0.65,
                         "ytick.major.width": 0.65, "savefig.facecolor": "white",
                         "figure.facecolor": "white", "lines.linewidth": 0.8})


def category_handles(context: Context) -> list[Patch]:
    """Make a shared legend whose threshold labels follow the actual configuration."""
    t = context.config["thresholds"]
    labels = {"core": f"Core (≥{100*t['core']:g}%)",
              "soft-core": f"Soft-core ({100*t['soft_core']:g}–<{100*t['core']:g}%)",
              "shell": f"Shell ({100*t['shell']:g}–<{100*t['soft_core']:g}%)",
              "cloud": f"Cloud (<{100*t['shell']:g}%)"}
    return [Patch(facecolor=COLORS[c], edgecolor="none", label=labels[c]) for c in CATEGORIES]


def new_figure(context: Context, height_mm: float, single: bool = False) -> plt.Figure:
    """Create an exact journal-column canvas; saving does not crop the physical width."""
    width = context.options["single_column_mm" if single else "double_column_mm"]
    return plt.figure(figsize=(width / 25.4, height_mm / 25.4), layout=None)


def save_figure(fig: plt.Figure, stem: str, context: Context) -> None:
    """Write vector and 600 dpi raster outputs with reproducible export metadata."""
    directory = context.outdir / "figures"
    directory.mkdir(exist_ok=True)
    for extension in context.options["formats"]:
        destination = directory / f"{stem}.{extension}"
        kwargs: dict = {"dpi": context.options["dpi"]}
        if extension == "pdf":
            kwargs["metadata"] = {"Creator": "pangenome_graph_figures.py", "CreationDate": None, "ModDate": None}
        elif extension == "svg":
            kwargs["metadata"] = {"Date": None, "Creator": "pangenome_graph_figures.py"}
        elif extension == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(destination, **kwargs)
        context.files.append(str(destination.relative_to(context.outdir)))
    LOG.info("Saved %s (%s)", stem, ", ".join(context.options["formats"]))
    plt.close(fig)


def export_source(panel: str, context: Context, node_ids: list[int], *,
                  positions: dict[int, np.ndarray] | None = None,
                  edges: list[tuple[int, int]] | None = None,
                  extras: list[dict] | None = None,
                  node_metadata: dict[int, dict] | None = None) -> None:
    """Export one CSV per panel, with separate graph and annotation evidence rows.

    Every node row joins to node_table.tsv by node_id. Edge endpoints use the
    same IDs. Aggregates are accompanied by the contributing individual nodes.
    """
    rows: list[dict] = []
    for n in sorted(set(node_ids)):
        row = context.table.loc[n]
        data = {key: row[key] for key in ("node_id", "size", "prevalence", "category", "degree", "paralog", "component_id")}
        data.update(record_type="node_measurement", evidence_type="graph_detected")
        if positions and n in positions:
            data.update(x=float(positions[n][0]), y=float(positions[n][1]))
        if node_metadata and n in node_metadata:
            data.update(node_metadata[n])
        rows.append(data)
        rows.append({"record_type": "annotation_label", "node_id": n,
                     "name": row["name"], "annotation": row["annotation"], "description": row["description"],
                     "panaroo_label": row["panaroo_label"], "evidence_type": "annotation_inferred"})
    for u, v in edges or []:
        rows.append({"record_type": "edge", "source_node_id": u, "target_node_id": v,
                     "edge_size": context.graph.edges[u, v]["size"],
                     "edge_support": context.graph.edges[u, v]["size"] / context.total,
                     "evidence_type": "graph_detected"})
    for extra in extras or []:
        rows.append({"evidence_type": "graph_detected", **extra})
    result = pd.DataFrame(rows)
    # Nullable IDs prevent float-formatted identifiers in mixed evidence tables.
    for key in ("node_id", "source_node_id", "target_node_id", "component_id"):
        if key in result:
            result[key] = pd.array(result[key], dtype="Int64")
    path = context.outdir / "figure_source_data" / f"{panel}.csv"
    result.to_csv(path, index=False)
    context.files.append(str(path.relative_to(context.outdir)))


def layout_graph(graph: nx.Graph, key: str, context: Context) -> dict[int, np.ndarray]:
    """Cache deterministic sparse layouts: igraph KK/DrL or Graphviz sfdp.

    The key includes input hash, exact topology, algorithm, software versions,
    and seed. No dense n-by-n graph matrices are created for large graphs.
    """
    ids = sorted(graph)
    if not ids:
        return {}
    if len(ids) == 1:
        return {ids[0]: np.zeros(2)}
    edge_list = sorted((min(u, v), max(u, v)) for u, v in graph.edges())
    requested = context.config["layout"]["algorithm"]
    try:
        import igraph as ig
        ig_version = ig.__version__
    except ImportError:
        ig = None
        ig_version = None
    if requested == "auto":
        algorithm = ("kk" if len(ids) <= 400 else "drl") if ig else "sfdp"
    else:
        algorithm = requested
    if algorithm in {"kk", "drl"} and ig is None:
        raise RuntimeError(f"Layout {algorithm} requires python-igraph; choose sfdp or auto")
    if algorithm == "kk" and len(ids) > 2000:
        raise ValueError("Kamada-Kawai is restricted to <=2000 nodes to avoid dense distance matrices; use DrL/sfdp")
    if algorithm not in {"kk", "drl", "sfdp"}:
        raise ValueError("layout.algorithm must be auto, kk, drl, or sfdp")
    if algorithm == "sfdp" and not shutil.which("sfdp"):
        raise RuntimeError("Install Graphviz (sfdp) or python-igraph; no large dense layout fallback is used")
    version = ig_version if algorithm != "sfdp" else subprocess.run(["sfdp", "-V"], capture_output=True, text=True).stderr.strip()
    identity = hashlib.sha256(json.dumps([CACHE_VERSION, context.fingerprint, ids, edge_list, algorithm, version, context.config["seed"]]).encode()).hexdigest()
    path = context.outdir / "cache" / f"layout_{key}_{identity[:16]}.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            xy = saved["xy"]
            if saved["node_ids"].tolist() != ids or not np.isfinite(xy).all():
                raise ValueError("Invalid layout cache")
        cache_hit = True
    else:
        LOG.info("Computing %s layout %s: %s nodes, %s edges", algorithm, key, len(ids), len(edge_list))
        rng = np.random.default_rng(context.config["seed"])
        if algorithm in {"kk", "drl"}:
            index = {n: i for i, n in enumerate(ids)}
            ig.set_random_number_generator(random.Random(context.config["seed"]))
            g = ig.Graph(n=len(ids), edges=[(index[u], index[v]) for u, v in edge_list], directed=False)
            initial = rng.normal(size=(len(ids), 2)).tolist()
            layout = g.layout_kamada_kawai(seed=initial) if algorithm == "kk" else g.layout_drl(seed=initial)
            xy = np.asarray(layout.coords, dtype=float)
        else:
            dot = "graph G {\n" + f'graph [start="{context.config["seed"]}", overlap="prism", pack="true"];\n'
            dot += 'node [shape="point", width="0.01", height="0.01", label=""];\n'
            dot += "\n".join(f"n{n};" for n in ids) + "\n"
            dot += "\n".join(f"n{u} -- n{v};" for u, v in edge_list) + "\n}\n"
            output = subprocess.run(["sfdp", "-Tplain"], input=dot, capture_output=True, text=True,
                                    timeout=context.config["layout"]["timeout_seconds"])
            if output.returncode:
                raise RuntimeError(f"Graphviz layout failed: {output.stderr[-1000:]}")
            coordinates = {int(parts[1][1:]): (float(parts[2]), float(parts[3]))
                           for line in output.stdout.splitlines() if (parts := line.split()) and parts[0] == "node"}
            xy = np.asarray([coordinates[n] for n in ids], dtype=float)
        xy -= xy.mean(axis=0)
        scale = max(float(np.ptp(xy, axis=0).max()), 1e-12)
        xy /= scale
        np.savez_compressed(path, node_ids=np.asarray(ids, dtype=np.int64), xy=xy)
        cache_hit = False
    context.layouts[key] = {"algorithm": algorithm, "software": version, "nodes": len(ids),
                            "edges": len(edge_list), "cache": str(path.relative_to(context.outdir)), "cache_hit": cache_hit,
                            "layout_edge_policy": "all edges in supplied subgraph; unweighted topology"}
    return dict(zip(ids, xy))


def gene_matches(gene_name: str, context: Context) -> list[int]:
    """Match case-insensitive exact name/annotation tokens; never phenotype descriptions."""
    target = gene_name.casefold()
    def matches(row: pd.Series) -> bool:
        terms = re.split(r"[;,]|~~~", f"{row['name']};{row['annotation']}")
        return any(term.strip().casefold() == target for term in terms)
    found = context.table[context.table.apply(matches, axis=1)].sort_values(["size", "node_id"], ascending=[False, True])
    return found.node_id.tolist()


def display_label(n: int, context: Context, maximum: int = 20) -> str:
    """Prefer an annotation gene symbol; retain the complete value in source data."""
    row = context.table.loc[n]
    alternatives = [term.strip() for term in str(row["annotation"]).split(";") if term.strip()]
    value = alternatives[0] if alternatives else str(row["name"])
    return value if len(value) <= maximum else value[:maximum - 1] + "…"


def repel_labels(texts: list, ax: plt.Axes, context: Context) -> None:
    """Use adjustText when available; deterministic display-space repulsion offline."""
    if not texts:
        return
    try:
        from adjustText import adjust_text
        # Iteration cap, rather than wall-clock stopping, preserves reproducibility.
        np.random.seed(context.config["seed"])
        adjust_text(texts, ax=ax, iter_lim=150, ensure_inside_axes=True,
                    arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.4})
        return
    except ImportError:
        context.warn("adjustText unavailable offline; using deterministic bounding-box label repulsion.")
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bounds = ax.get_window_extent(renderer)
    occupied = []
    for text in texts:
        origin = np.asarray(ax.transData.transform(text.get_position()))
        chosen = origin.copy()
        found = False
        for attempt in range(300):
            if attempt:
                angle = attempt * 2.399963229728653
                radius = 3.5 * np.sqrt(attempt)
                chosen = origin + radius * np.array([np.cos(angle), np.sin(angle)])
            chosen = np.clip(chosen, [bounds.x0 + 35, bounds.y0 + 8], [bounds.x1 - 35, bounds.y1 - 8])
            text.set_position(ax.transData.inverted().transform(chosen))
            bbox = text.get_window_extent(renderer).expanded(1.07, 1.16)
            if not any(bbox.overlaps(other) for other in occupied):
                occupied.append(bbox)
                found = True
                break
        if not found:
            text.set_visible(False)
            context.warn("A label was suppressed to avoid overlap; all annotations remain in source data.")
        elif np.linalg.norm(chosen - origin) > 4:
            destination = ax.transData.inverted().transform(chosen)
            ax.annotate("", xy=ax.transData.inverted().transform(origin), xytext=destination,
                        arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.4}, zorder=4)


def select_edges(graph: nx.Graph, context: Context, *, cap: bool = False) -> list[tuple[int, int]]:
    """Filter display edges by cohort support without changing measured topology."""
    cutoff = context.options["min_edge_support"]
    edges = [(min(u, v), max(u, v)) for u, v, a in graph.edges(data=True) if a["size"] / context.total >= cutoff]
    edges.sort(key=lambda e: (-graph.edges[e]["size"], e[0], e[1]))
    if cap and len(edges) > context.options["max_overview_edges"]:
        context.warn(f"Overview edges capped deterministically at the {context.options['max_overview_edges']} highest-support edges.")
        edges = edges[:context.options["max_overview_edges"]]
    return edges


def node_areas(ids: list[int], context: Context, small: bool) -> np.ndarray:
    """Bound marker areas with a monotone log1p(prevalence) transform."""
    prevalence = context.table.loc[ids, "prevalence"].to_numpy(dtype=float)
    return (0.7 + 2.0 * np.log1p(prevalence) / np.log(2)) if small else (8 + 28 * np.log1p(prevalence) / np.log(2))


def draw_graph(ax: plt.Axes, graph: nx.Graph, positions: dict, context: Context, *,
               edges: list[tuple[int, int]] | None = None, label_ids: list[int] | None = None,
               small: bool = False, centers: list[int] | None = None, rasterized: bool = False,
               label_overrides: dict[int, str] | None = None) -> None:
    """Draw sparse edge segments and category layers, with edge support as width."""
    if not graph:
        ax.text(0.5, 0.5, "No nodes meet these criteria", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return
    edges = list(graph.edges()) if edges is None else edges
    if edges:
        segments = np.asarray([[positions[u], positions[v]] for u, v in edges])
        supports = np.asarray([graph.edges[u, v]["size"] / context.total for u, v in edges])
        lines = LineCollection(segments, colors="#53585D", linewidths=(0.12 + 0.6*supports) if small else (0.25 + 1.4*supports),
                               alpha=0.20 if small else 0.36, rasterized=rasterized, zorder=1)
        ax.add_collection(lines)
    # Cloud below the more prevalent clusters, always in fixed order.
    for category in reversed(CATEGORIES):
        ids = [n for n in sorted(graph) if context.table.loc[n, "category"] == category]
        if ids:
            xy = np.asarray([positions[n] for n in ids])
            ax.scatter(xy[:, 0], xy[:, 1], s=node_areas(ids, context, small), c=COLORS[category],
                       linewidths=0.1 if not small else 0, edgecolors="white" if not small else "none",
                       alpha=0.9, rasterized=rasterized, zorder=2)
    if centers:
        xy = np.asarray([positions[n] for n in centers if n in positions])
        if len(xy):
            ax.scatter(xy[:, 0], xy[:, 1], s=70, facecolors="none", edgecolors="#202328", linewidths=1.0, zorder=3)
    xy = np.asarray(list(positions.values()))
    span = np.ptp(xy, axis=0)
    margin = np.maximum(span * 0.10, 0.04)
    ax.set_xlim(xy[:, 0].min()-margin[0], xy[:, 0].max()+margin[0], auto=True)
    ax.set_ylim(xy[:, 1].min()-margin[1], xy[:, 1].max()+margin[1], auto=True)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_axis_off()
    texts = [ax.text(*positions[n], (label_overrides or {}).get(n, display_label(n, context)), fontsize=7, zorder=5,
                     bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none", "pad": 0.6})
             for n in label_ids or [] if n in positions]
    repel_labels(texts, ax, context)


def plot_composition(context: Context, single: bool = False) -> plt.Figure:
    """Figure A: log-frequency histogram and four-category cluster counts."""
    fig = new_figure(context, 130 if single else 80, single=single)
    if single:
        grid = fig.add_gridspec(2, 1, left=.22, right=.96, top=.87, bottom=.24, hspace=.75)
    else:
        grid = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], left=0.09, right=0.98, top=0.83, bottom=0.28, wspace=0.36)
    ax = fig.add_subplot(grid[0]); counts_ax = fig.add_subplot(grid[1])
    draw_composition(ax, context)
    counts = context.table.category.value_counts().reindex(CATEGORIES, fill_value=0)
    y = np.arange(4)
    counts_ax.barh(y, counts.values, color=[COLORS[c] for c in CATEGORIES], height=0.62)
    counts_ax.set_yticks(y, [c.title() for c in CATEGORIES]); counts_ax.invert_yaxis()
    counts_ax.set_xlim(0, max(counts.max()*1.37, 1)); counts_ax.set_xlabel("Gene clusters")
    for i, count in enumerate(counts):
        counts_ax.text(count + counts.max()*0.025, i, f"{count:,}", va="center", fontsize=7.5)
    counts_ax.set_title("Pangenome composition", loc="left")
    fig.suptitle(f"{context.config['dataset_label']} · {context.total:,} genomes", fontsize=9 if single else 11, y=0.98)
    fig.legend(handles=category_handles(context), loc="lower center", ncol=1 if single else 2, frameon=False, bbox_to_anchor=(0.5, 0.015 if single else 0.025))
    edges = np.linspace(0, 1, context.options["histogram_bins"]+1)
    extras = [{"record_type": "category_count", "category": c, "gene_count": int(counts[c])} for c in CATEGORIES]
    for c in CATEGORIES:
        hist, _ = np.histogram(context.table.loc[context.table.category == c, "prevalence"], bins=edges)
        extras.extend({"record_type": "frequency_bin", "category": c, "bin_left": edges[i], "bin_right": edges[i+1], "gene_count": int(count)} for i, count in enumerate(hist))
    if not single:
        export_source("panel_a_composition", context, context.table.node_id.tolist(), extras=extras)
    return fig


def draw_composition(ax: plt.Axes, context: Context, compact: bool = False) -> None:
    """Draw the frequency histogram with thresholds and inclusive final bin."""
    bins = np.linspace(0, 1, context.options["histogram_bins"]+1)
    samples = [context.table.loc[context.table.category == c, "prevalence"].to_numpy() for c in CATEGORIES]
    ax.hist(samples, bins=bins, stacked=True, color=[COLORS[c] for c in CATEGORIES], edgecolor="none", log=True)
    ax.set_xlim(0, 1); ax.set_ylim(bottom=0.8)
    # Plain numerals avoid log-tick superscripts below the 7 pt minimum.
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.set_xticks([0, .25, .5, .75, 1], ["0", "25", "50", "75", "100"])
    ax.set_xlabel("Genomes carrying cluster (%)"); ax.set_ylabel("Gene clusters (log scale)")
    for cutoff in context.config["thresholds"].values():
        ax.axvline(cutoff, color="#666666", linestyle=(0, (2, 2)), lw=0.5, zorder=0)
    ax.set_title("Gene frequency", loc="left")
    if compact:
        counts = context.table.category.value_counts().reindex(CATEGORIES, fill_value=0)
        inset = ax.inset_axes([0.26, 0.55, 0.45, 0.37])
        inset.barh(np.arange(4), counts, color=[COLORS[c] for c in CATEGORIES], height=.65)
        inset.invert_yaxis(); inset.set_xlim(0, max(counts.max()*1.8, 1))
        inset.set_xticks([]); inset.set_yticks([])
        for i, (c, count) in enumerate(counts.items()):
            inset.text(count + counts.max()*.04, i, f"{count:,}", va="center", fontsize=7)
        for spine in inset.spines.values():
            spine.set_visible(False)


def backbone_data(context: Context) -> tuple[nx.Graph, dict, list, list]:
    """Select core/soft-core nodes, retaining all induced edges for layout."""
    ids = context.table.loc[context.table.category.isin(["core", "soft-core"]), "node_id"].tolist()
    graph = context.graph.subgraph(ids)
    positions = layout_graph(graph, "backbone", context)
    edges = select_edges(graph, context)
    ranked = sorted(ids, key=lambda n: (-context.graph.degree(n), n))[:context.options["top_labels"]]
    requested = list(itertools.chain.from_iterable(gene_matches(gene, context) for gene in context.options["label_genes"]))
    labels = list(dict.fromkeys([n for n in requested if n in graph] + ranked))
    return graph, positions, edges, labels


def plot_backbone(context: Context) -> plt.Figure:
    """Figure B: core/soft-core induced graph with support-weighted display edges."""
    graph, positions, edges, labels = backbone_data(context)
    fig = new_figure(context, 118); ax = fig.add_axes([0.035, .15, .93, .70])
    draw_graph(ax, graph, positions, context, edges=edges, label_ids=labels, small=True)
    fig.suptitle("Core and soft-core backbone", fontsize=11, y=.97)
    fig.text(.5, .90, f"{len(graph):,} clusters · {len(edges):,}/{graph.number_of_edges():,} edges shown · support ≥{context.options['min_edge_support']:.0%}", ha="center", fontsize=8)
    fig.legend(handles=category_handles(context)[:2], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5,.045))
    fig.text(.5, .025, "Node area: log1p(prevalence). Edge width: genome support. Labels: annotated as named.", ha="center", fontsize=7)
    export_source("panel_b_backbone", context, list(graph), positions=positions, edges=edges,
                  node_metadata={n: {"node_area_pt2": float(a), "label_requested": n in labels} for n, a in zip(sorted(graph), node_areas(sorted(graph), context, True))},
                  extras=[{"record_type": "display_policy", "min_edge_support": context.options["min_edge_support"], "edges_total": graph.number_of_edges(), "edges_shown": len(edges)}])
    return fig


def overview_data(context: Context) -> tuple[dict, list]:
    """Keep every node and lay out all topology, filtering only rendered edges."""
    return layout_graph(context.graph, "overview", context), select_edges(context.graph, context, cap=True)


def draw_overview(ax: plt.Axes, positions: dict, edges: list, context: Context) -> None:
    """Use rasterized sparse artists above 20k nodes; editable text stays vector."""
    raster = len(context.graph) > context.options["rasterize_above_nodes"]
    # Support filtering is the deliberate scalable alternative to datashader.
    # No nodes are sampled, and the complete graph determines the cached layout.
    draw_graph(ax, context.graph, positions, context, edges=edges, small=True, rasterized=raster)
    context.renderers["overview"] = "support-filtered Matplotlib; rasterized node/edge artists" if raster else "support-filtered Matplotlib vector artists"
    ax.text(.02, .015, f"All {len(context.graph):,} nodes\n{len(edges):,}/{context.graph.number_of_edges():,} edges shown", transform=ax.transAxes, fontsize=7,
            va="bottom", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9, "pad": 2})


def plot_overview(context: Context) -> plt.Figure:
    """Figure C: complete-node graph overview with explicit rendering policy."""
    positions, edges = overview_data(context)
    fig = new_figure(context, 125); ax = fig.add_axes([.025, .17, .95, .72])
    draw_overview(ax, positions, edges, context)
    fig.suptitle("Full pangenome graph", fontsize=11, y=.97)
    fig.legend(handles=category_handles(context), loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5,.055))
    fig.text(.5, .018, f"Display edges: support ≥{context.options['min_edge_support']:.0%}. Layout coordinates are arbitrary; no genomic distance scale.", ha="center", fontsize=7)
    export_source("panel_c_overview", context, list(context.graph), positions=positions, edges=edges,
                  extras=[{"record_type": "display_policy", "min_edge_support": context.options["min_edge_support"], "edges_total": context.graph.number_of_edges(), "edges_shown": len(edges), "all_nodes_retained": True}])
    return fig


def hotspot_table(context: Context) -> tuple[pd.DataFrame, float]:
    """Rank unfiltered core degrees; accessory adjacency is a junction proxy."""
    table = context.table.loc[context.table.category == "core"].copy()
    if table.empty:
        return table, float("nan")
    median = float(table.degree.median())
    mad = float(np.median(np.abs(table.degree - median)))
    # A one-degree minimum scale is explicit because MAD is often zero for
    # a mostly linear backbone. This is a descriptive cutoff, not a p-value.
    cutoff = median + context.options["hotspot_mad_multiplier"] * max(1.0, 1.4826 * mad)
    table["accessory_neighbours"] = [sum(context.table.loc[v, "category"] in {"shell", "cloud"} for v in context.graph[n]) for n in table.node_id]
    table["candidate_junction"] = (table.degree > cutoff) & (table.accessory_neighbours > 0)
    table = table.sort_values(["degree", "accessory_neighbours", "node_id"], ascending=[False, False, True])
    table["rank"] = np.arange(1, len(table)+1)
    table["degree_cutoff"] = cutoff
    return table, cutoff


def zoom_graph(centers: list[int], radius: int, key: str, context: Context) -> tuple[nx.Graph, dict, list[int], dict]:
    """Build a radius-limited ego union, with a deterministic connected display cap."""
    distances: dict[int, int] = {}
    for center in centers:
        for n, d in nx.single_source_shortest_path_length(context.graph, center, cutoff=radius).items():
            distances[n] = min(d, distances.get(n, d))
    ranked = sorted(distances, key=lambda n: (n not in centers, distances[n], -context.table.loc[n, "size"], -context.graph.degree(n), n))
    cap = max(context.options["max_zoom_nodes"], len(centers))
    chosen = ranked[:cap]
    graph = context.graph.subgraph(chosen)
    positions = layout_graph(graph, key, context)
    labels = list(dict.fromkeys(centers + [n for n in chosen if context.table.loc[n, "annotation"]]))[:context.options["max_zoom_labels"]]
    meta = {"centers": centers, "radius": radius, "ego_nodes_total": len(distances), "nodes_displayed": len(chosen),
            "nodes_omitted": len(distances)-len(chosen), "edges_displayed": graph.number_of_edges(),
            "node_policy": "centers; then hop distance, decreasing prevalence, decreasing degree, node_id",
            "edge_policy": "all induced edges between displayed nodes; no support filter"}
    if meta["nodes_omitted"]:
        context.warn(f"Zoom {key}: showing {len(chosen)}/{len(distances)} nodes at radius {radius}; deterministic cap recorded in source data.")
    context.zooms[key] = meta
    return graph, positions, labels, meta


def draw_hotspot_ranking(ax: plt.Axes, table: pd.DataFrame, cutoff: float, context: Context,
                         compact: bool = False) -> pd.DataFrame:
    """Draw the lollipop ranking without interpreting gene names as phenotypes."""
    top = table.head(6 if compact else context.options["hotspot_top_n"])
    if top.empty:
        ax.text(.5,.5,"No core clusters",ha="center",transform=ax.transAxes)
        return top
    y = np.arange(len(top))
    ax.hlines(y, 0, top.degree, color="#CBD2D7", lw=1)
    ax.scatter(top.degree, y, color=COLORS["core"], s=18, zorder=3)
    ax.set_yticks(y, [display_label(n, context, 15 if compact else 24) for n in top.node_id])
    ax.invert_yaxis(); ax.set_xlabel("Degree in unfiltered graph")
    ax.set_xlim(0, max(top.degree.max()*1.22, cutoff*1.12))
    ax.axvline(cutoff, color="#777777", linestyle="--", lw=.7)
    ax.set_title("Core junction candidates", loc="left")
    for yi, (_, row) in zip(y, top.iterrows()):
        ax.text(row.degree+top.degree.max()*.02, yi, str(row.degree), va="center", fontsize=7)
    return top


def plot_hotspots(context: Context) -> plt.Figure:
    """Figure D: ranked core degrees and top-three radius-two neighbourhoods."""
    table, cutoff = hotspot_table(context)
    table.to_csv(context.outdir / "hotspot_candidates.tsv", sep="\t", index=False)
    fig = new_figure(context, 154)
    ax = fig.add_axes([.17,.51,.76,.32])
    top = draw_hotspot_ranking(ax, table, cutoff, context)
    fig.suptitle("Accessory-adjacency junction candidates", fontsize=11, y=.97)
    fig.text(.5,.91,"Graph-derived ranking; high degree does not establish an insertion event",ha="center",fontsize=8)
    ids = set(table.node_id)
    extras = [{"record_type": "hotspot_rank", "node_id": int(row.node_id), "rank": int(row['rank']),
               "accessory_neighbours": int(row.accessory_neighbours), "candidate_junction": bool(row.candidate_junction),
               "degree_cutoff": cutoff} for _, row in table.iterrows()]
    selected = table[table.candidate_junction].head(3)
    if len(selected) < 3:
        selected = table.head(3)
    all_edges = set()
    for index, n in enumerate(selected.node_id):
        graph, positions, labels, meta = zoom_graph([n], 2, f"hotspot_{n}", context)
        zoom_ax = fig.add_axes([.015 + index*.33, .11, .31, .28])
        draw_graph(zoom_ax, graph, positions, context, label_ids=labels[:7], centers=[n])
        zoom_ax.set_title(f"Annotated as {display_label(n, context, 15)}\n{len(graph)}/{meta['ego_nodes_total']} nodes · radius 2",fontsize=7.5, pad=3)
        ids.update(graph); all_edges.update((min(u,v),max(u,v)) for u,v in graph.edges())
        extras.append({"record_type": "zoom_policy", "zoom_id": f"hotspot_{n}", **meta})
        extras.extend({"record_type": "zoom_node_position", "zoom_id": f"hotspot_{n}", "node_id": v,
                       "x": float(positions[v][0]), "y": float(positions[v][1]), "is_center": v == n} for v in graph)
    fig.legend(handles=category_handles(context), loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5,.025))
    export_source("panel_d_hotspots", context, sorted(ids), edges=sorted(all_edges), extras=extras)
    return fig


def neighbourhood_data(gene_name: str, radius: int, context: Context) -> tuple[nx.Graph, dict, list, dict] | None:
    """Resolve all exact annotation matches and cache their ego-union layout."""
    matches = gene_matches(gene_name, context)
    if not matches:
        context.warn(f"Focus gene {gene_name!r} has no exact name/annotation match; this does not demonstrate biological absence.")
        return None
    if len(matches) > 1:
        context.warn(f"Focus gene {gene_name!r} matches {len(matches)} clusters; all matches are retained as centres.")
    slug = re.sub(r"[^A-Za-z0-9_-]", "_", gene_name)
    return zoom_graph(matches, radius, f"focus_{slug}", context)


def plot_neighbourhood(gene_name: str, radius: int = 2, *, context: Context) -> plt.Figure | None:
    """Figure E: ego neighbourhood of a gene name (annotation match) or gene cluster.

    Exact matches in `name` and semicolon-separated `annotation` are searched.
    All matches are included, with ambiguity and display caps made explicit.
    An unmatched gene produces a warning rather than an invented plot.
    """
    data = neighbourhood_data(gene_name, radius, context)
    if data is None:
        return None
    graph, positions, labels, meta = data
    fig = new_figure(context, 112); ax = fig.add_axes([.025,.17,.95,.66])
    draw_graph(ax, graph, positions, context, label_ids=labels, centers=meta["centers"],
               label_overrides={n:gene_name for n in meta["centers"]})
    fig.suptitle(f"Gene neighbourhood · annotated as {gene_name}",fontsize=11,y=.97)
    fig.text(.5,.89,f"Radius {radius} · {len(meta['centers'])} annotation match(es) · {len(graph)}/{meta['ego_nodes_total']} nodes displayed",ha="center",fontsize=8)
    fig.legend(handles=category_handles(context),loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(.5,.055))
    fig.text(.5,.018,"Outlined nodes: annotation matches. All displayed induced edges are shown. No phenotypic inference.",ha="center",fontsize=7)
    slug = re.sub(r"[^A-Za-z0-9_-]", "_", gene_name)
    export_source(f"panel_e_neighbourhood_{slug}",context,list(graph),positions=positions,edges=list(graph.edges()),
                  node_metadata={n:{"is_center":n in meta["centers"], "displayed_label":gene_name if n in meta["centers"] else display_label(n,context),
                                    "label_requested":n in labels} for n in graph},
                  extras=[{"record_type":"zoom_policy", "focus_gene":gene_name, **meta}])
    return fig


def component_table(context: Context) -> pd.DataFrame:
    """Summarize connected components of the complete, unfiltered graph."""
    rows=[]
    for cid, frame in context.table.groupby("component_id",sort=True):
        counts=frame.category.value_counts()
        rows.append({"component_id":int(cid), "node_count":len(frame),
                     "edge_count":context.graph.subgraph(frame.node_id).number_of_edges(),
                     **{f"{c.replace('-','_')}_count":int(counts.get(c,0)) for c in CATEGORIES},
                     "accessory_only":bool(frame.category.isin(["shell","cloud"]).all()),
                     "evidence_type":"graph_detected"})
    return pd.DataFrame(rows)


def draw_component_distribution(ax: plt.Axes, components: pd.DataFrame) -> None:
    """Draw exact discrete component sizes on logarithmic axes."""
    sizes=components.node_count.value_counts().sort_index()
    ax.scatter(sizes.index, sizes.values, s=18, color="#53585D", zorder=3)
    ax.vlines(sizes.index, .75, sizes.values, color="#B5BDC4", lw=.7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(bottom=.75)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.set_xlabel("Clusters per connected component"); ax.set_ylabel("Components (log scale)")
    ax.set_title("Connected-component sizes",loc="left")


def plot_components(context: Context) -> plt.Figure:
    """Figure F: complete-graph components and three largest accessory-only examples."""
    components=component_table(context)
    components.to_csv(context.outdir / "component_table.tsv",sep="\t",index=False)
    fig=new_figure(context,142)
    ax=fig.add_axes([.12,.57,.8,.27]); draw_component_distribution(ax,components)
    fig.suptitle("Pangenome connected components",fontsize=11,y=.97)
    fig.text(.5,.91,f"{len(components):,} unfiltered components · {int(components.accessory_only.sum()):,} accessory-only",ha="center",fontsize=8)
    selected=components[components.accessory_only].sort_values(["node_count","component_id"],ascending=[False,True]).head(3)
    extras=[{"record_type":"component_summary", **row} for row in components.to_dict("records")]
    shown_edges=set()
    for index,row in enumerate(selected.itertuples()):
        ids=context.table.loc[context.table.component_id==row.component_id,"node_id"].tolist()
        graph=context.graph.subgraph(ids)
        positions=layout_graph(graph,f"component_{row.component_id}",context)
        sub=fig.add_axes([.025+index*.33,.15,.30,.24])
        draw_graph(sub,graph,positions,context)
        sub.set_title(f"Component {row.component_id}\n{len(graph)} accessory clusters",fontsize=7.5)
        shown_edges.update((min(u,v),max(u,v)) for u,v in graph.edges())
        extras.extend({"record_type":"component_node_position", "node_id":n,"component_id":row.component_id,
                       "x":float(positions[n][0]),"y":float(positions[n][1])} for n in graph)
    if selected.empty:
        fig.text(.5,.31,"No complete-graph component contains only shell/cloud clusters.",ha="center",fontsize=8)
    fig.legend(handles=category_handles(context),loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(.5,.045))
    fig.text(.5,.015,"Disconnected accessory components are graph observations; mobile-element or island identity is unverified.",ha="center",fontsize=7)
    export_source("panel_f_components",context,list(context.graph),edges=sorted(shown_edges),extras=extras)
    return fig


def panel_letter(ax: plt.Axes, letter: str, offset: tuple[float,float] = (-.16,1.08)) -> None:
    """Add consistent bold lowercase letters at final print size."""
    ax.text(*offset,letter,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top",clip_on=False)


def plot_multipanel(context: Context, focus_genes: list[str]) -> plt.Figure:
    """Compose A–F with native artists, a shared legend, and exact double-column width."""
    fig=new_figure(context,232)
    grid=fig.add_gridspec(3,2,left=.105,right=.97,bottom=.14,top=.92,wspace=.43,hspace=.60)
    axa=fig.add_subplot(grid[0,0]); draw_composition(axa,context,compact=True); panel_letter(axa,"a")
    axb=fig.add_subplot(grid[0,1]); graph,pos,edges,labels=backbone_data(context)
    draw_graph(axb,graph,pos,context,edges=edges,label_ids=labels[:3],small=True)
    axb.set_title("Core / soft-core backbone",loc="left"); panel_letter(axb,"b")
    axc=fig.add_subplot(grid[1,0]); positions,edges=overview_data(context); draw_overview(axc,positions,edges,context)
    axc.set_title("Full pangenome",loc="left"); panel_letter(axc,"c")
    axd=fig.add_subplot(grid[1,1]); hot,cutoff=hotspot_table(context)
    draw_hotspot_ranking(axd,hot,cutoff,context,compact=True); panel_letter(axd,"d")
    axe=fig.add_subplot(grid[2,0]); focus=None
    for gene in focus_genes:
        if gene_matches(gene,context):
            focus=gene; break
    if focus:
        graph,pos,labels,meta=neighbourhood_data(focus,context.options["neighbourhood_radius"],context)
        draw_graph(axe,graph,pos,context,label_ids=labels[:5],centers=meta["centers"],
                   label_overrides={n:focus for n in meta["centers"]})
        axe.set_title(f"Annotated as {focus} · radius {meta['radius']}",loc="left",fontsize=9)
    else:
        axe.set_axis_off(); axe.text(.5,.5,"No requested annotation match",ha="center",va="center",transform=axe.transAxes)
        axe.set_title("Gene neighbourhood",loc="left")
    panel_letter(axe,"e")
    axf=fig.add_subplot(grid[2,1])
    if context.options["component_analysis"]:
        draw_component_distribution(axf,component_table(context))
    else:
        axf.set_axis_off(); axf.text(.5,.5,"Component analysis disabled",ha="center",transform=axf.transAxes)
    panel_letter(axf,"f")
    fig.suptitle(f"{context.config['organism']} · {context.config['dataset_label'].lower()}\n{context.total:,} genomes · {len(context.graph):,} gene clusters",fontsize=11,y=.985)
    fig.legend(handles=category_handles(context),loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(.5,.065),columnspacing=2)
    fig.text(.5,.037,f"Graph geometry and category counts are detected from the graph. Gene labels are annotation-inferred.\nOverview/backbone edge display: support ≥{context.options['min_edge_support']:.0%}; layouts use all induced edges. No phenotypic claims.",ha="center",fontsize=7)
    index=pd.DataFrame([{"composite_panel":letter,"source_csv":source,"evidence_type":"graph_detected"}
                        for letter,source in [("a","panel_a_composition.csv"),("b","panel_b_backbone.csv"),("c","panel_c_overview.csv"),
                                               ("d","panel_d_hotspots.csv"),("e",f"panel_e_neighbourhood_{re.sub(r'[^A-Za-z0-9_-]','_',focus)}.csv" if focus else "no match"),
                                               ("f","panel_f_components.csv" if context.options["component_analysis"] else "disabled")]])
    index.to_csv(context.outdir/"figure_source_data"/"multipanel_source_index.csv",index=False)
    return fig


def software_versions() -> dict[str, str | None]:
    """Record exact runtime versions, including missing optional accelerators."""
    result: dict[str,str|None] = {"python":sys.version.split()[0]}
    for name in ["networkx","igraph","matplotlib","pandas","numpy","datashader","adjustText","PyYAML","scipy","seaborn","tqdm"]:
        try:
            result[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name]=None
    result["graphviz"] = subprocess.run(["sfdp","-V"],capture_output=True,text=True).stderr.strip() if shutil.which("sfdp") else None
    return result


def write_legends_and_methods(context: Context, denominator_source: str, versions: dict, focus_genes: list[str]) -> None:
    """Generate manuscript-ready methods from the actual run, not guessed versions."""
    t=context.config["thresholds"]; o=context.options
    counts=context.table.category.value_counts().reindex(CATEGORIES,fill_value=0)
    algorithms=sorted({info['algorithm'] for info in context.layouts.values()})
    version_text=", ".join(f"{name} {versions[name]}" for name in ["python","networkx","matplotlib","pandas","numpy","scipy","PyYAML"] if versions[name])
    methods=(f"The Panaroo GML graph was analysed using {version_text}. The supplied graph contained {len(context.graph):,} gene-cluster nodes and "
             f"{context.graph.number_of_edges():,} undirected adjacency edges. The cohort denominator ({context.total:,} genomes) was obtained from {denominator_source}. "
             f"Gene-cluster prevalence was defined as the node size attribute divided by the cohort size. Clusters were classified as core (≥{100*t['core']:g}%), "
             f"soft-core (≥{100*t['soft_core']:g}% and <{100*t['core']:g}%), shell (≥{100*t['shell']:g}% and <{100*t['soft_core']:g}%), or cloud (<{100*t['shell']:g}%). "
             f"Counts were {counts['core']:,}, {counts['soft-core']:,}, {counts['shell']:,}, and {counts['cloud']:,}, respectively. "
             f"Node sizes and edge supports were checked against membership cardinalities. Graphs were positioned using {', '.join(algorithms)} layouts "
             f"({versions.get('graphviz') if 'sfdp' in algorithms else 'igraph '+str(versions.get('igraph'))}) with a fixed seed of {context.config['seed']}; coordinates were cached. "
             f"Layouts used all edges in each supplied graph/subgraph and unweighted topology. Backbone and overview displays retained edges supported by ≥{o['min_edge_support']:.0%} "
             f"of the cohort; measurements of degree and connected components used the complete graph. Overview node/edge layers were rasterized when exceeding "
             f"{o['rasterize_above_nodes']:,} nodes. Core junction candidates were defined by degree greater than the core median plus "
             f"{o['hotspot_mad_multiplier']:g} × max(1, 1.4826 × median absolute deviation), together with at least one shell/cloud neighbour. "
             f"This descriptive rule has no inferential significance threshold. Neighbourhoods comprised radius-{o['neighbourhood_radius']} ego unions around exact "
             f"name/annotation matches; displays were capped at {o['max_zoom_nodes']} nodes by hop distance, prevalence, degree and node ID. "
             "Display omissions and exact edge lists were exported as source data. Gene symbols and descriptions were treated as annotation-inferred labels, "
             "without asserting experimentally established resistance, virulence, mobile-element identity, or insertion events. The original Panaroo version, "
             "annotation program/version, and upstream QC settings were not encoded reliably in the supplied files and must be supplied from the original analysis records. "
             f"Figures used an Okabe-Ito palette, ≥7 pt sans-serif typography, editable PDF/SVG text, and {o['dpi']} dpi raster exports.")
    (context.outdir/"methods.md").write_text("# Methods for this run\n\n"+methods+"\n")
    legends=f"""# Figure legends and evidence conventions

**a. Gene frequency and composition.** Prevalence is graph node size / {context.total:,} genomes. Histogram categories are stacked on a logarithmic y-axis. Counts refer to gene clusters, not gene copies or genomes. Boundaries are inclusive at {100*t['core']:g}%, {100*t['soft_core']:g}%, and {100*t['shell']:g}% for core, soft-core, and shell, respectively.

**b. Core/soft-core backbone.** All qualifying nodes are included. Edge width represents size / cohort size; displayed edges meet the {o['min_edge_support']:.0%} support cutoff. Node area is a bounded affine transform of log1p(prevalence). Labels identify clusters annotated as the stated name; labels are selected by full-graph degree or explicit requests. Layout is based on all induced edges.

**c. Full graph overview.** Every cluster is retained. The layout uses the complete graph; rendered edges meet the {o['min_edge_support']:.0%} support cutoff and are capped, if necessary, at {o['max_overview_edges']:,} highest-support edges with deterministic ID tie-breaking. Rasterized graph layers preserve editable legend/text. Coordinates are arbitrary and convey no genomic or phylogenetic distance. Isolated plotted nodes may have incident edges omitted from display.

**d. Core junction candidates.** Lollipops show full-graph degrees of ranked core clusters. The dashed line is the median-plus-MAD descriptive cutoff specified in Methods; candidate status also requires accessory adjacency. Top-three candidate ego graphs include radius-two neighbours; all induced edges between displayed nodes are shown. Each zoom reports displayed/total ego nodes. Degree does not establish an accessory insertion site; paralogy, assembly boundaries and graph cleaning may affect adjacency. Cluster names are annotation-inferred.

**e. Gene neighbourhood.** Requested genes: {', '.join(focus_genes) or 'none'}. Nodes outlined in black are exact name/annotation matches, annotated as the requested gene. Multiple matches are included and reported. All induced edges between displayed nodes are shown regardless of global support cutoff. No match is an annotation-search result and does not establish biological absence. Category colors are graph-derived. No resistance or virulence phenotype is inferred.

**f. Components.** Component sizes are computed from the complete, unfiltered graph. Standalone examples show the largest complete-graph components made entirely of shell/cloud nodes, rather than components obtained by removing core nodes. Disconnected accessory components could warrant follow-up, but their mobile-element or genomic-island identity has not been demonstrated.

**Combined figure.** Panels a-f summarize their standalone counterparts. The combined junction panel shows the top six ranks; three detailed hotspot zooms are available in standalone d. Panel e uses the first requested gene with an exact match. Panel f shows the component-size distribution; example accessory-only components are available in standalone f. All underlying values are in the panel CSV files indexed by `multipanel_source_index.csv`.

**Source-data convention.** `node_measurement`, `edge`, rank, component, and layout rows have `evidence_type=graph_detected`. `annotation_label` rows have `evidence_type=annotation_inferred`. Join `node_id` (or edge endpoint IDs) to `node_table.tsv`; annotation-derived fields in that master table carry `label_evidence_type=annotation_inferred`. Layout and display selections are computational summaries of graph data, not experimental results. `panaroo_label` preserves the input GML label in addition to the GML integer ID. CSV record types are heterogeneous by design; blank cells are not zero values.
"""
    (context.outdir/"figure_legends.md").write_text(legends)


def validate_source_data(context: Context) -> dict[str, int]:
    """Verify plotted IDs/edge endpoints and evidence labels against the master table."""
    valid=set(context.table.node_id)
    result={}
    for path in sorted((context.outdir/"figure_source_data").glob("panel_*.csv")):
        data=pd.read_csv(path, low_memory=False)
        for column in ["node_id","source_node_id","target_node_id"]:
            if column in data:
                values=set(data[column].dropna().astype(int))
                if not values <= valid:
                    raise AssertionError(f"Source IDs absent from node_table.tsv: {path.name}, {column}")
        if not set(data.evidence_type) <= {"graph_detected","annotation_inferred"}:
            raise AssertionError(f"Invalid evidence_type in {path}")
        result[path.name]=len(data)
    return result


def run_pipeline(args: argparse.Namespace) -> Context:
    """Inventory, inspect, categorize, plot, validate traceability, and write provenance."""
    inputs=detect_inputs(Path.cwd())
    print("Auto-detected inputs in current directory:\n"+json.dumps(inputs,indent=2),flush=True)
    path=choose_gml(inputs,args.gml)
    print(f"Selected graph: {path}",flush=True)
    config=read_config(args.config)
    if config["dataset_label"] is None:
        cohort="Clinical" if "clinical" in path.stem.lower() else "Environmental" if "env" in path.stem.lower() else "Pangenome"
        config["dataset_label"]=f"{cohort} pangenome" if cohort != "Pangenome" else cohort
    if args.min_edge_support is not None:
        if not 0 <= args.min_edge_support <= 1:
            raise ValueError("--min-edge-support must be a fraction in [0,1]")
        config["figures"]["min_edge_support"]=args.min_edge_support
    focus=[value.strip() for value in args.focus_genes.split(",") if value.strip()] if args.focus_genes is not None else config["figures"]["focus_genes"]
    config["figures"]["focus_genes"]=focus
    outdir=args.outdir.resolve(); outdir.mkdir(parents=True,exist_ok=True)
    for directory in ["cache","figure_source_data","figures"]:
        (outdir/directory).mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(),logging.FileHandler(outdir/"pipeline.log",mode="w")],force=True)
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    started=time.monotonic()
    fingerprint=sha256_file(path)
    graph,inspection=load_gml(path,outdir/"cache",config,fingerprint)
    if not graph:
        raise ValueError("Graph has no nodes")
    total,denominator_source=infer_total_genomes(inspection,config["total_genomes"])
    # Inspection is printed and persisted before category or figure calculation.
    write_inspection(graph,inspection,total,denominator_source,path,fingerprint,outdir/"00_gml_inspection.txt")
    table=make_node_table(graph,total,config["thresholds"])
    table.to_csv(outdir/"node_table.tsv",sep="\t",index=False)
    context=Context(graph,table,total,config,outdir,fingerprint)
    if "union" in denominator_source:
        context.warn("Genome count inferred from observed membership union; genomes with no retained genes cannot be counted. Verify total_genomes against upstream records.")
    for kind,expected,keys,missing in [("node",EXPECTED_NODE_KEYS,inspection.node_keys,inspection.node_missing),
                                      ("edge",EXPECTED_EDGE_KEYS,inspection.edge_keys,inspection.edge_missing)]:
        absent=expected-keys
        if absent:
            context.warn(f"Expected {kind} attributes absent: {', '.join(sorted(absent))}")
        if missing:
            context.warn(f"Some {kind} attributes are missing on individual records: {dict(missing)}")
    if inspection.node_size_disagreements or inspection.edge_size_disagreements:
        raise ValueError("size/member-cardinality discrepancies prevent validated prevalence/support; see inspection report")
    if any(a["size"] > total for _,_,a in graph.edges(data=True)):
        raise ValueError("Edge support exceeds the genome denominator")
    if inspection.derived_node_sizes or inspection.derived_edge_sizes:
        context.warn(f"Missing size values derived from membership: {inspection.derived_node_sizes} nodes, {inspection.derived_edge_sizes} edges.")
    if "fallback" in inspection.probe:
        context.warn(inspection.probe)
    versions=software_versions()
    if versions["igraph"] is None:
        context.warn("python-igraph unavailable offline; auto layouts use installed Graphviz sfdp.")
    if versions["datashader"] is None and len(graph)>config["figures"]["rasterize_above_nodes"]:
        context.warn("Datashader unavailable offline; full overview uses support-filtered sparse edges and rasterized Matplotlib node/edge layers.")
    configure_style(context)
    if config["statistics"]["enabled"]:
        from pangenome_statistics import run_statistics
        context.statistics = run_statistics(context, inspection.isolate_names)
        context.table.to_csv(outdir/"node_table.tsv", sep="\t", index=False)
    (outdir/"config_resolved.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    (outdir/"software_versions.json").write_text(json.dumps(versions,indent=2)+"\n")
    with (outdir/"requirements_runtime.txt").open("w") as handle:
        for name,version in versions.items():
            if version and name not in {"python","graphviz"}:
                handle.write(f"{name}=={version}\n")
    functions=[("fig_a_composition",plot_composition),("fig_b_backbone",plot_backbone),
               ("fig_c_overview",plot_overview),("fig_d_hotspots",plot_hotspots)]
    for stem,function in tqdm(functions,desc="Figures A-D"):
        save_figure(function(context),stem,context)
    save_figure(plot_composition(context,single=True),"fig_a_composition_single_column",context)
    match_rows=[]
    for gene in tqdm(focus,desc="Gene neighbourhoods"):
        matches=gene_matches(gene,context)
        match_rows.extend({"requested_gene":gene,"node_id":n,"name":table.loc[n,"name"],"annotation":table.loc[n,"annotation"],"evidence_type":"annotation_inferred"} for n in matches)
        figure=plot_neighbourhood(gene,config["figures"]["neighbourhood_radius"],context=context)
        if figure is not None:
            save_figure(figure,"fig_e_neighbourhood_"+re.sub(r"[^A-Za-z0-9_-]","_",gene),context)
    pd.DataFrame(match_rows,columns=["requested_gene","node_id","name","annotation","evidence_type"]).to_csv(outdir/"focus_gene_matches.tsv",sep="\t",index=False)
    if config["figures"]["component_analysis"]:
        save_figure(plot_components(context),"fig_f_components",context)
    save_figure(plot_multipanel(context,focus),"pangenome_multipanel",context)
    source_counts=validate_source_data(context)
    write_legends_and_methods(context,denominator_source,versions,focus)
    if context.statistics:
        from pangenome_statistics import write_statistical_methods
        write_statistical_methods(context)
    manifest={"input_gml":str(path),"input_sha256":fingerprint,"input_inventory":inputs,
              "pipeline_sha256":sha256_file(Path(__file__)),"total_genomes":total,"denominator_source":denominator_source,
              "nodes":len(graph),"edges":graph.number_of_edges(),"connected_components":nx.number_connected_components(graph),
              "category_counts":table.category.value_counts().reindex(CATEGORIES,fill_value=0).astype(int).to_dict(),
              "config":config,"software_versions":versions,"layouts":context.layouts,"zooms":context.zooms,
              "statistics":context.statistics,
              "statistics_module_sha256":sha256_file(Path(__file__).with_name("pangenome_statistics.py")),
              "renderers":context.renderers,"figure_files":context.files,"source_csv_rows":source_counts,
              "warnings":context.warnings,"elapsed_seconds":round(time.monotonic()-started,2)}
    (outdir/"run_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    summary=[f"Generated results in {outdir}",f"Genomes: {total:,}; nodes: {len(graph):,}; edges: {graph.number_of_edges():,}; components: {nx.number_connected_components(graph):,}",
             "Categories: "+", ".join(f"{c}={manifest['category_counts'][c]:,}" for c in CATEGORIES),
             f"Figure exports: {len([f for f in context.files if f.startswith('figures/')])}; source CSVs: {len(source_counts)}",
             "Master table: node_table.tsv; inspection: 00_gml_inspection.txt",
             "Methods, legends, exact versions, resolved config, and run manifest written.",
             "Warnings:"]+[f"- {warning}" for warning in context.warnings]
    if context.statistics:
        summary.insert(3, f"Louvain communities: {context.statistics['communities']:,}; tested pairs: {context.statistics['tested_pairs']:,}; BH-significant positive pairs: {context.statistics['significant_pairs']:,}")
    (outdir/"summary.txt").write_text("\n".join(summary)+"\n")
    print("\n"+"\n".join(summary),flush=True)
    return context


def main() -> int:
    """Single CLI entry point for the complete visualization pipeline."""
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gml",type=Path,help="GML input; defaults to canonical or uniquely named clinical graph")
    parser.add_argument("--outdir",type=Path,default=Path("results"),help="Output directory (default: results)")
    parser.add_argument("--config",type=Path,default=Path("config.yaml") if Path("config.yaml").exists() else None,help="YAML configuration")
    parser.add_argument("--focus-genes",help="Comma-separated exact gene names; overrides configuration")
    parser.add_argument("--min-edge-support",type=float,help="Display cutoff: edge size / total genomes (0-1)")
    args=parser.parse_args()
    try:
        run_pipeline(args)
        return 0
    except (ValueError,FileNotFoundError,RuntimeError,subprocess.TimeoutExpired) as exc:
        LOG.error("Pipeline stopped: %s",exc)
        print(f"ERROR: {exc}",file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
