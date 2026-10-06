#!/usr/bin/env python3
from common import build_parser, context_from_args, supplementary_tables

if __name__ == "__main__":
    args = build_parser("Generate supplementary tables and dictionaries.").parse_args()
    supplementary_tables(context_from_args(args))
