#!/bin/bash
# Count variants in PLINK2 .snplist file

SNPLIST=$1

if [ -z "$SNPLIST" ]; then
    echo "Error: No file specified"
    echo "Usage: $0 <snplist_file>"
    exit 1
fi

if [ ! -f "$SNPLIST" ]; then
    echo "Error: File not found: $SNPLIST"
    exit 1
fi

echo "Counting variants in: $SNPLIST"

VARIANT_COUNT=$(wc -l < "$SNPLIST")

echo "Total variants: $VARIANT_COUNT"