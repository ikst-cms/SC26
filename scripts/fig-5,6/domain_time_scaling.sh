#!/bin/bash

for i in 1 2 6
do
	mkdir -p quaternary-node-$i
	chmod +x test-node-$i.sh
	cp test-node-$i.sh quaternary-node-$i/job.sh
	cp -r quaternary-scratch quaternary-node-$i/
	cd quaternary-node-$i
	qsub -V job.sh
	cd ../
done


output_csv="$(pwd)/mc_summary.csv"
input_file="MC.log"

# Write header if CSV doesn't exist
if [ ! -f "$output_csv" ]; then
    echo "node,model,supercell_size,avg_step_secs,total_metro_steps,total_metro_time_hrs,total_flips,time_per_flip_secs,avg_nbrs,p_length,p_sites" > "$output_csv"
fi

for node_dir in quaternary-node-*; do
    [ -d "$node_dir" ] || continue
    for core_dir in "$node_dir"/quaternary-core-*; do
        [ -d "$core_dir" ] || continue
        model=$(basename "$core_dir")
        for size in 100 1000 10000 100000; do
            size_dir="$core_dir/$size"
            [ -d "$size_dir" ] || continue
            file="$size_dir/$input_file"
            [ -f "$file" ] || continue

            avg_nbrs=$(grep "Average number of neighbors" "$file" | awk -F": " '{print $2}' | head -1)
            p_length=$(grep "Perturbation length uesd in calculating energy difference" "$file" | awk -F": " '{print $2}' | head -1)
            p_sites=$(grep "Approximate Number of sites within perturbed region" "$file" | awk -F": " '{print $2}' | awk '{print $1}' | head -1)

            # Get all Metropolis times in hours and average per step, output in seconds
            nlines=$(grep -c "Time taken for 10 Metropolis equilibriation steps" "$file")
            sumhrs=$(grep "Time taken for 10 Metropolis equilibriation steps" "$file" | awk '{sum+=$(NF-1)} END {print sum}')
            if [ "$nlines" -gt 0 ]; then
                avg_step_secs=$(awk "BEGIN {print ($sumhrs / ($nlines * 10)) * 3600}")
            else
                avg_step_secs=""
            fi

            # Extract total metropolis steps (sum across all lines)
            total_metro_steps=$(grep "Metropolis equilibriation steps" "$file" \
                                | grep -oE '[0-9]+ Metropolis' \
                                | grep -oE '[0-9]+' \
                                | awk '{sum += $1} END {if (NR > 0) print sum; else print "N/A"}')

            # Extract total metropolis time in hours (sum across all lines)
            total_metro_time_hrs=$(grep "Metropolis equilibriation steps" "$file" \
                                   | grep -oE '[0-9]+\.[0-9]+ hrs' \
                                   | grep -oE '[0-9]+\.[0-9]+' \
                                   | awk '{sum += $1} END {if (NR > 0) printf "%.6f\n", sum; else print "N/A"}')

            # Calculate total_flips from custom_list_size file
            custom_file="$size_dir/custom_list_size"
            if [ -f "$custom_file" ]; then
                total_flips=$(awk 'NF {for(i=1;i<=NF;i++) sum+=$i; has_data=1} END {if(has_data) print sum; else print "N/A"}' "$custom_file")
            else
                total_flips="N/A"
            fi

            # Calculate time per flip in seconds
            if [ "$total_flips" != "N/A" ] && [ -n "$total_metro_time_hrs" ] && [ "$total_metro_time_hrs" != "N/A" ]; then
                time_per_flip_secs=$(awk "BEGIN {print ($total_metro_time_hrs * 3600) / $total_flips}")
            else
                time_per_flip_secs="N/A"
            fi

            # Append extracted data to CSV
            echo "$(basename "$node_dir"),$model,$size,$avg_step_secs,$total_metro_steps,$total_metro_time_hrs,$total_flips,$time_per_flip_secs,$avg_nbrs,$p_length,$p_sites" >> "$output_csv"
            echo "Data extracted for $node_dir/$model/$size and appended."
        done
    done
done
echo "Data extraction complete. Summary saved to $output_csv."