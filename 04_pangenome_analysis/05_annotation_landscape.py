#!/usr/bin/env python3
from common import build_parser, context_from_args, annotation_landscape

if __name__ == "__main__":
    args = build_parser("Summarize functional annotation and length landscape.").parse_args()
    annotation_landscape(context_from_args(args))
