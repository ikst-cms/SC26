#!/bin/bash

output_csv="mc_scaling_summary.csv"

# CSV Header - folder_name removed
echo "atoms,cores,metropolis_steps,metropolis_time_hrs" > "$output_csv"

# Loop through all log files (1.log, 2.log, 4.log, etc.)
for log_file in [0-9]*.log; do
  # Skip if no matching files
  if [ ! -f "$log_file" ]; then
    continue
  fi
  
  # Extract core count from filename (e.g., 256.log -> 256)
  cores=$(echo "$log_file" | grep -oE '^[0-9]+')
  
  # Extract number of atoms (active sites)
  atoms=$(grep "active sites out of" "$log_file" | head -n 1 | awk '{print $2}')
  
  # Extract and SUM all Metropolis steps
  metro_steps=$(grep "Metropolis equilibriation steps" "$log_file" \
                | grep -oE '[0-9]+ Metropolis' \
                | grep -oE '[0-9]+' \
                | awk '{sum += $1} END {if (NR > 0) print sum}')
  
  # Extract and SUM all Metropolis times
  metro_time=$(grep "Metropolis equilibriation steps" "$log_file" \
               | grep -oE '[0-9]+\.[0-9]+ hrs' \
               | grep -oE '[0-9]+\.[0-9]+' \
               | awk '{sum += $1} END {if (NR > 0) printf "%.6f\n", sum}')

  # Append all extracted values to the CSV (no folder_name)
  echo "${atoms:-N/A},${cores:-N/A},${metro_steps:-N/A},${metro_time:-N/A}" >> "$output_csv"
done

echo "Extraction complete. Results saved to $output_csv"