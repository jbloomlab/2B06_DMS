#!/bin/bash

set -e

RESULTSDIR="results"
mkdir -p $RESULTSDIR


declare -a nbs=(
                "analysis_notebook.ipynb"
                "neut_analysis.ipynb"
                )


for nb in "${nbs[@]}"
do
    echo "Running $nb"

    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $nb

    jupyter nbconvert \
        --output-dir $RESULTSDIR \
        --to markdown \
        $nb
done
