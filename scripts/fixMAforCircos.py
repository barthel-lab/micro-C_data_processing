#!/bin/python
import sys
import os

# Get input file from the command line argument
input_file = sys.argv[1]

# Generate output file name by changing the extension to .fix.txt
output_file = os.path.splitext(input_file)[0] + ".itinMA.txt"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    pairs = []
    for line in infile:
        fields = line.strip().split("\t")
        pairs.append(fields)

    # Ensure we process pairs of lines
    for i in range(0, len(pairs), 2):
        if i + 1 < len(pairs):  # Avoid index error if an odd number of lines
            link_name = f"link{i//2 + 1}"
            # For the first chromosome
            chr1, start1, end1 = pairs[i][1], pairs[i][2], pairs[i][3]
            outfile.write(f"{link_name}\t{chr1}\t{start1}\n")
            # For the second chromosome
            chr2, start2, end2 = pairs[i+1][1], pairs[i+1][2], pairs[i+1][3]
            outfile.write(f"{link_name}\t{chr2}\t{start2}\t{end2}\n")

print(f"Processed data saved to {output_file}")

