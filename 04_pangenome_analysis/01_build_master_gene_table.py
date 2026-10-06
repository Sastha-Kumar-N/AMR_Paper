#!/usr/bin/env python3
from common import build_parser, context_from_args, build_master_gene_table

if __name__ == "__main__":
    args = build_parser("Build the master gene-cluster table.").parse_args()
    build_master_gene_table(context_from_args(args))
