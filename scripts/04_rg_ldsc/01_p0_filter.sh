#!/bin/bash

# Hard-coded file paths to summary statistics (.ma) files
# Header of .ma file: SNP	A1	A2	Freq	b	se	p	N
input_file="/.../SUMSTAT_INPUT.ma"
output_file="/.../SUMSTAT_OUT.ma"

# Recommended threshold for p-values (adjust as needed) https://github.com/bulik/ldsc/issues/144
threshold=1e-300

# Use awk to filter out rows where the p-value (column 7) is below the threshold.
# Assumes the file is tab-delimited and the header is on the first line.
awk -v thresh="$threshold" 'BEGIN { FS="\t"; OFS="\t" }
NR==1 { print; next }
{
    # Convert the p-value to a number
    p_val = $7 + 0
    if (p_val >= thresh) print
}' "$input_file" > "$output_file"

echo "Rows with a p-value below $threshold have been removed."
echo "Filtered output saved to: $output_file"
