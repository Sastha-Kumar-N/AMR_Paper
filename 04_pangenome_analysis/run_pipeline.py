#!/usr/bin/env python3
from common import build_parser, context_from_args, run_all

if __name__ == "__main__":
    args = build_parser("Run the complete whole-pangenome publication workflow.").parse_args()
    run_all(context_from_args(args), args)
