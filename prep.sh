#!/bin/bash

#nohup ./prep.sh GBR > GBR.log 2>&1 &
# ./prep.sh CEU 4
POP=$1

DIR="/home/ailina/DAIseg.grch37.real"

process_chromosome() {
    i=$1

    json="${DIR}/${POP}.YRI.jsons/${POP}.YRI.grch37.chr${i}.json"

    [[ -f "$json" ]] && \
    python daiseg.py restrict_1kG -json "$json" -threads 16 && \
    python daiseg.py callability -json "$json" -threads 16 && \
    python daiseg.py main.prep -json "$json" -threads 16 && \
    { 
        if [[ ! -f "$json" ]]; then
            echo "ERROR: JSON file $json not found"
        else
            echo "ERROR: Processing failed for chr$i"
        fi
        return 1
    }

    return 0
}

export -f process_chromosome
export POP
export JOBS

chromosomes=({22..1})


for chr in "${chromosomes[@]}"; do
    echo " Starting chromosome $chr at $(date)"
    process_chromosome "$chr" || {
        echo "ERROR: Failed on chromosome $chr"
        exit 1
    }

    echo " Completed chromosome $chr at $(date)"
done

