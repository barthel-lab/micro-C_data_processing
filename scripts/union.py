#!/bin/python
import sys

def load_names(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip().split()[0][1:] for line in f)

# Ensure correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python union.py <input_file1> <input_file2> <output_file>")
    sys.exit(1)

# Load file paths from command line arguments
file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = sys.argv[3]

# Load read names from both files
names1 = load_names(file1)
names2 = load_names(file2)

# Find union
all_names = names1.union(names2)

# Write union to output file
with open(output_file, 'w') as f:
    for name in sorted(all_names):
        f.write(name + "\n")

print(f"Union saved to {output_file}")