#!/usr/bin/env python3
from common import build_parser, context_from_args, accessory_genome_structure

if __name__ == "__main__":
    args = build_parser("Analyze accessory-genome ordination and clustering.").parse_args()
    accessory_genome_structure(context_from_args(args))
