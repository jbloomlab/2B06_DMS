# Mutational antigenic profiling 2B06 antibody

## Overview

This directory contains the analysis of the mutational antigenic profiling of A/WSN/1933 (H1N1) influenza virus hemagglutinin escape from the neutralizing antibody d045-051310-2B06, from here on refered to as 2B06.

Antibody 2B06 was provided by Patrick Wilson, and first characterized by [Dunand et al. 2015](https://www.jci.org/articles/view/74374). 

Selection of virus libraries with 2B06 was performed by Lauren Gentles as was subsequent analysis.

## Results
The easiest way to look at the results is to view the Markdown rendering the Jupyter notebook at [results/analysis_notebook.md](results/analysis_notebook.md).

## Running the analysis
The analysis is contained in the Python Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
To run the notebook and generate the Markdown summary at [results/analysis_notebook.md](results/analysis_notebook.md), use the bash script [run_notebook.bash](run_notebook.bash) by running:

    ./run_notebook.bash

To submit and run this on the Hutch cluster, you do:

    sbatch -c 4 run_notebook.bash

# Input data
Input data are in the [./data/](data) subdirectory.

Input files include the [WSN reference sequence](./data/WSN-HA.fasta) and [samplelist.csv](./data/samplelist.csv) which provides the R1 location as well as relavent information about each sample.
