#!/usr/bin/env python3
from common import build_parser, context_from_args, per_genome_gene_content

if __name__ == "__main__":
    args = build_parser("Calculate per-genome gene-content summaries.").parse_args()
    per_genome_gene_content(context_from_args(args))
