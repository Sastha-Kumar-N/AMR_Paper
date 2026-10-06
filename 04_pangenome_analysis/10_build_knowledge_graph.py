#!/usr/bin/env python3
from common import build_parser, context_from_args, build_knowledge_graph

if __name__ == "__main__":
    args = build_parser("Build the pangenome knowledge graph.").parse_args()
    build_knowledge_graph(context_from_args(args))
