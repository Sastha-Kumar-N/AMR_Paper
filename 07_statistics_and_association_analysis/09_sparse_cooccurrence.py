#!/usr/bin/env python3
from common import build_parser, context_from_args, sparse_cooccurrence

if __name__ == "__main__":
    args = build_parser("Run sparse filtered gene co-occurrence analysis.").parse_args()
    sparse_cooccurrence(
        context_from_args(args),
        max_genes=args.cooccurrence_max_genes,
        min_prevalence=args.cooccurrence_min_prevalence,
        max_prev_frac=args.cooccurrence_max_prevalence_frac,
        min_shared=args.cooccurrence_min_shared,
    )
