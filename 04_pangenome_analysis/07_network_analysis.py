#!/usr/bin/env python3
from common import build_parser, context_from_args, network_analysis

if __name__ == "__main__":
    args = build_parser("Analyze Panaroo gene-neighborhood graph topology.").parse_args()
    network_analysis(context_from_args(args))
