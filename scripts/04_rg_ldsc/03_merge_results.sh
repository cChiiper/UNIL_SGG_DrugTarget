#!/bin/bash

# Directory containing log files
DIRECTORY="/.../"

# Output file
OUTPUT_FILE="/.../pw_genetic_correlation_results.tsv"

# Remove existing output file if it exists
if [[ -f "$OUTPUT_FILE" ]]; then
    rm "$OUTPUT_FILE"
fi

# Write the header (18 columns: 12 from the summary table + gcov, gcov_se + 4 new heritability columns)
echo -e "p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se\tgcov\tgcov_se\th2_p1\th2_se_p1\th2_p2\th2_se_p2" > "$OUTPUT_FILE"

# Process each log file in the directory
for file in "$DIRECTORY"/*.log; do

  # Extract genetic covariance values (supports scientific notation)
  gcov=$(grep "Total Observed scale gencov:" "$file" | head -n1 | \
         sed -E 's/.*gencov:[[:space:]]+(-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?)[[:space:]]+\(.*/\1/')
  gcov_se=$(grep "Total Observed scale gencov:" "$file" | head -n1 | \
            sed -E 's/.*\((-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?)\).*/\1/')

  # Extract heritability values for phenotype 1
  # Using -A 2 so that we capture the line with h2 (2 lines after the header)
  h2_line_p1=$(grep -A 2 "Heritability of phenotype 1" "$file" | grep "Total Observed scale h2:" | head -n1)
  h2_p1=$(echo "$h2_line_p1" | sed -E 's/.*h2:[[:space:]]+(-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?).*/\1/')
  h2_se_p1=$(echo "$h2_line_p1" | sed -E 's/.*\((-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?)\).*/\1/')

  # Extract heritability values for phenotype 2 (using the header with "2/2")
  h2_line_p2=$(grep -A 2 "Heritability of phenotype 2/2" "$file" | grep "Total Observed scale h2:" | head -n1)
  h2_p2=$(echo "$h2_line_p2" | sed -E 's/.*h2:[[:space:]]+(-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?).*/\1/')
  h2_se_p2=$(echo "$h2_line_p2" | sed -E 's/.*\((-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?)\).*/\1/')

  # Extract the genetic correlation summary table from after "Summary of Genetic Correlation Results"
  # until "Analysis finished".  Skip the header row of the table.
  awk '/Summary of Genetic Correlation Results/{flag=1; next} /Analysis finished/{flag=0} flag' "$file" |
  tail -n +2 |
  awk -v gcov_val="$gcov" -v gcov_se_val="$gcov_se" \
      -v h2_p1_val="$h2_p1" -v h2_se_p1_val="$h2_se_p1" \
      -v h2_p2_val="$h2_p2" -v h2_se_p2_val="$h2_se_p2" \
      -v OFS='\t' '{
    # Clean up whitespace from the line
    sub(/^[ \t]+/, "", $0); sub(/[ \t]+$/, "", $0);
    if ($0 != "") {
      # Extract only the filename for p1 and p2 by splitting on "/"
      split($1, p1_arr, "/");
      split($2, p2_arr, "/");
      # Print all 12 columns from the summary table followed by gcov, gcov_se,
      # h2_p1, h2_se_p1, h2_p2, and h2_se_p2
      print p1_arr[length(p1_arr)], p2_arr[length(p2_arr)], $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, gcov_val, gcov_se_val, h2_p1_val, h2_se_p1_val, h2_p2_val, h2_se_p2_val
    }
  }' >> "$OUTPUT_FILE"

done

echo "Data extraction complete. Results saved to $OUTPUT_FILE."
