# Mutational antigenic profiling 2B06 antibody

## Overview

This directory contains the analysis of the mutational antigenic profiling of A/WSN/1933 (H1N1) influenza virus hemagglutinin escape from the neutralizing antibody d045-051310-2B06, from here on refered to as 2B06.

Antibody 2B06 was provided by Patrick Wilson, and first characterized by [Dunand et al. 2015](https://www.jci.org/articles/view/74374). 

Selection of virus libraries with 2B06 was performed by Lauren Gentles as was subsequent analysis.

## Results
The easiest way to look at the results is to view the Markdown rendering the Jupyter notebook at [results/analysis_notebook.md](results/analysis_notebook.md).

In addition, the results can be visualized on the structure by opening [map_on_struct.ipynb](map_on_struct.ipynb) on [mybinder](https://mybinder.org/) by entering this link: https://mybinder.org/v2/gh/jbloomlab/2B06_DMS/master?urlpath=%2Fapps%2Fmap_on_struct.ipynb

## Running the analysis
The analysis is contained in the Python Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
To run the notebook and generate the Markdown summary at [results/analysis_notebook.md](results/analysis_notebook.md), use the bash script [run_nbs.bash](run_nbs.bash) by running:

    ./run_nbs.bash

To submit and run this on the Hutch cluster, you do:

    sbatch -c 4 run_nbs.bash

# Input data
Input data are in the [./data/](data) subdirectory.

- [./data/WSN-HA.fasta](./data/WSN-HA.fasta) contains the wildtype sequence of the A/WSN/1933 H1N1 hemagglutinin used in the experimental virus library.

- [./data/samplelist.csv](./data/samplelist.csv) provides the R1 location as well as relavent information about each sample.

- [./data/H1toH3renumber.csv](./data/H1toH3_renumber.csv) is a CSV file that mapes the numbering from the sequential numbering of the WSN HA protein to the commonly used H3 nubering scheme. The sequential numbering is in the *original* column, and the H3 number is in the *new* column.

- [./data/H3_site_to_PDB_1rvx.csv](./data/H3_site_to_PDB_1rvx.csv) is a file that maps the sites in WSN HA H3 numbering to the equivalent sites in the PR8 HA PDB structure [1rvx](https://www.rcsb.org/structure/1RVX). This file has columns site, which is the H3 site, the PDB chain, which is the chain in PDB file [1rvx](https://www.rcsb.org/structure/1RVX), and the PDB site, which is the site in the PDB.

- [./data/neut_data/2B06_neut_analysis.xlsx](./data/neut_data/2B06_neut_analysis.xlsx) is a file containing the raw data from the GFP based neutralization assay measuring the neutralization potency of 2B06 against wildtype WSN GFP virus.

- [./data/neut_data/2B06_neut_022019.yaml](./data/neut_data/2B06_neut_022019.yaml) is a file that provides the location of the neutralization data file and the name of the samples to [neut_analysis.ipynb](neut_analysis.ipynb) which generates the neutralization curves from this data.
