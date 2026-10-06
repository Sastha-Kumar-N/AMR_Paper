#!/usr/bin/env python3
from common import build_parser, context_from_args, pangenome_composition

if __name__ == "__main__":
    args = build_parser("Calculate pangenome composition and occupancy summaries.").parse_args()
    pangenome_composition(context_from_args(args))
