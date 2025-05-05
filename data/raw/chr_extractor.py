#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.fasta"
    exit 1
fi

input="$1"

# Use awk to split the FASTA file
awk '
    # When a header line is encountered (starts with ">")
    /^>/ {
        # Close the previous file if it exists
        if (outfile) close(outfile)
        # Remove the leading ">" and take the first token as the identifier
        split($0, tokens, " ")
        header = substr(tokens[1], 2)
        outfile = header ".fasta"  # e.g. Chr1.fasta, Chr2.fasta, etc.
    }
    # Write the current line to the current output file
    { print >> outfile }
' "$input"
