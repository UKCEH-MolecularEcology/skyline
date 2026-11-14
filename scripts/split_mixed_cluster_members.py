#!/usr/bin/env python3

import csv
import sys

# very large clusters (e.g. 90% identity) have member lists longer than csv's default field limit (131072 chars)
_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_limit)
        break
    except OverflowError:
        _limit //= 10

# Get input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Open the input CSV for reading and the output CSV for writing
with open(input_file, newline="") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.DictReader(infile)  # Read input as a dictionary per row
    # Add a new field 'Member_Type' to the existing fieldnames
    fieldnames = reader.fieldnames + ["Member_Type"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()  # Write the new header to output

    # Iterate through each row in the input
    for row in reader:
        # Split the 'Members' field by semicolon
        members = row["Members"].split(";")

        # For each member, generate a new row
        for member in members:
            member = member.strip()
            new_row = row.copy()  # Copy the existing row
            new_row["Members"] = member  # Replace the 'Members' field with a single member

            # Identify whether the member is short read (ASV) or long read (contig)
            if member.startswith("ASV"):
                new_row["Member_Type"] = "amplicon"
            else:
                new_row["Member_Type"] = "long_read"

            writer.writerow(new_row)  # Write the updated row
