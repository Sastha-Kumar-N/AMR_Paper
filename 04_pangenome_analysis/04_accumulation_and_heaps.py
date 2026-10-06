#!/usr/bin/env python3
from common import build_parser, context_from_args, accumulation_and_heaps

if __name__ == "__main__":
    args = build_parser("Run pangenome accumulation and Heaps' law analyses.").parse_args()
    accumulation_and_heaps(context_from_args(args), bootstrap=args.accumulation_bootstrap)
