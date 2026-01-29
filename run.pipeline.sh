#!/bin/bash

# Set the number of parallel jobs (adjust based on your system)
JOBS=6

# Create a function to process each chromosome
process_chromosome() {
    i=$1
    echo "=== Starting chromosome $i ==="
    
    json="./CHB.YRI.jsons/CHB.YRI.grch37.chr${i}.json"
    
    # Check if JSON file exists
    if [[ ! -f "$json" ]]; then
        echo "ERROR: JSON file $json not found"
        return 1
    fi
    
    echo "Processing $json"
    
    # Step 1: restrict_1kG
    echo "  Step 1: restrict_1kG"
    python daiseg.py restrict_1kG -json "$json" -threads 4
    if [[ $? -ne 0 ]]; then
        echo "ERROR: restrict_1kG failed for chr$i"
        return 1
    fi
    
    # Step 2: callability
    echo "  Step 2: callability"
    python daiseg.py callability -json "$json" -threads 4
    if [[ $? -ne 0 ]]; then
        echo "ERROR: callability failed for chr$i"
        return 1
    fi
    
    # Step 3: main.prep
    echo "  Step 3: main.prep"
    python daiseg.py main.prep -json "$json" -threads 4
    if [[ $? -ne 0 ]]; then
        echo "ERROR: main.prep failed for chr$i"
        return 1
    fi
    
    echo "=== Finished chromosome $i successfully ==="
    return 0
}

# Export function to make it available to parallel
export -f process_chromosome

# Create array of chromosomes (22 to 1)
chromosomes=({22..1})

echo "Starting parallel processing of ${#chromosomes[@]} chromosomes with $JOBS parallel jobs"
echo "Chromosomes to process: ${chromosomes[@]}"

# Run with GNU Parallel
parallel --jobs $JOBS \
         --progress \
         --joblog daiseg_joblog.txt \
         --results daiseg_results \
         --tag \
         process_chromosome ::: "${chromosomes[@]}"

# Check exit status
if [[ $? -eq 0 ]]; then
    echo ""
    echo "========================================"
    echo "All chromosomes processed successfully!"
    echo "========================================"
    
    # Show summary
    echo ""
    echo "Job log summary:"
    cat daiseg_joblog.txt | tail -n +2 | awk '{print "Chromosome", $1, ":", ($7==0) ? "SUCCESS" : "FAILED", "Exit code:", $7}'
else
    echo ""
    echo "========================================"
    echo "Some jobs failed! Check the output above."
    echo "========================================"
fi
