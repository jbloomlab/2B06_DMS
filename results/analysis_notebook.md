<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Mutational-antigenic-profiling-of-2B06-antibody-selection." data-toc-modified-id="Mutational-antigenic-profiling-of-2B06-antibody-selection.-1">Mutational antigenic profiling of 2B06 antibody selection.</a></span><ul class="toc-item"><li><span><a href="#Overview" data-toc-modified-id="Overview-1.1">Overview</a></span></li><li><span><a href="#Import-Python-packages" data-toc-modified-id="Import-Python-packages-1.2">Import Python packages</a></span></li><li><span><a href="#Process-deep-sequencing-data" data-toc-modified-id="Process-deep-sequencing-data-1.3">Process deep sequencing data</a></span><ul class="toc-item"><li><span><a href="#Get-information-about-samples" data-toc-modified-id="Get-information-about-samples-1.3.1">Get information about samples</a></span></li><li><span><a href="#Run-dms2_batch_bcsubamp" data-toc-modified-id="Run-dms2_batch_bcsubamp-1.3.2">Run <code>dms2_batch_bcsubamp</code></a></span></li><li><span><a href="#Look-at-summary-plots" data-toc-modified-id="Look-at-summary-plots-1.3.3">Look at summary plots</a></span></li></ul></li><li><span><a href="#Renumber-codon-counts-to-H3-numbering" data-toc-modified-id="Renumber-codon-counts-to-H3-numbering-1.4">Renumber codon counts to H3 numbering</a></span></li><li><span><a href="#Compute-differential-selection" data-toc-modified-id="Compute-differential-selection-1.5">Compute differential selection</a></span></li><li><span><a href="#Logo-plots-of-differential-selection" data-toc-modified-id="Logo-plots-of-differential-selection-1.6">Logo plots of differential selection</a></span></li></ul></li></ul></div>
# Mutational antigenic profiling of 2B06 antibody selection.

## Overview
The goal of this analysis notebook is to determine the comprehensive set of single amino-acid mutations in influenza A/WSN/33 hemagglutinin that allow for virus escape from the antibody D045-051310-2B06.
This antibody was provided by Patrick Wilson and previously characterized as a broadly cross-reactive HA [antibody](https://www.ncbi.nlm.nih.gov/pubmed/25689254). Deep mutational scanning experiments were performed by Lauren Gentles _in_ _vitro_. Selections were performed on MDCK cells and percent infectivity after antibody selection was determined by qPCR.

## Import Python packages


```python
import os
import pandas
from IPython.display import display, HTML
import dms_tools2
from dms_tools2.ipython_utils import showPDF
import glob

print("Using dms_tools2 version {0}".format(dms_tools2.__version__))

# results will go in this directory
resultsdir = './results/' 
if not os.path.isdir(resultsdir):
    os.mkdir(resultsdir)
    
# CPUs to use, should not exceed the number you request with slurm
ncpus = 14

# do we use existing results or generate everything new?
use_existing = 'no'

fastqdir = os.path.join(resultsdir, 'FASTQ_files/')
if not os.path.isdir(fastqdir):
    os.mkdir(fastqdir)
    
renumberedcountsdir = os.path.join(resultsdir, 'renumberedcounts')
if not os.path.isdir(renumberedcountsdir):
    os.mkdir(renumberedcountsdir)
```

    Using dms_tools2 version 2.4.12


## Process deep sequencing data

### Get information about samples
Here we define the samples by reading in a *.csv* file containing the selection type, library, percent infectivity, and R1 path. We then combind the selection type and the library to create the sample name. For example, the WSN library 1 selected with antibody 2B06 gets the name, '2B06-Lib.1'. We also include a mock unselected sample for comparison to identified enriched mutations in the 2B06 selected sample. The WSN-plasmid encodes the wildtype HA sequence and will served as an error control in this analysis. There are two separate error control samples in this analysis since two different sequencing runs were required to produce the data.


```python
samples = pandas.read_csv('data/samplelist.csv')
samples['name'] = samples.selection + '-' + samples.library

display(HTML(samples.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>percent_infectivity</th>
      <th>R1</th>
      <th>name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2B06</td>
      <td>Lib.1</td>
      <td>0.17</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
      <td>2B06-Lib.1</td>
    </tr>
    <tr>
      <td>Mock</td>
      <td>Lib.1</td>
      <td></td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
      <td>Mock-Lib.1</td>
    </tr>
    <tr>
      <td>WSN-plasmid</td>
      <td></td>
      <td></td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
      <td>WSN-plasmid-</td>
    </tr>
    <tr>
      <td>2B06</td>
      <td>Lib.2</td>
      <td>0.10</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
      <td>2B06-Lib.2</td>
    </tr>
    <tr>
      <td>Mock</td>
      <td>Lib.2</td>
      <td></td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
      <td>Mock-Lib.2</td>
    </tr>
    <tr>
      <td>2B06</td>
      <td>Lib.3</td>
      <td>0.14</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
      <td>2B06-Lib.3</td>
    </tr>
    <tr>
      <td>Mock</td>
      <td>Lib.3</td>
      <td></td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
      <td>Mock-Lib.3</td>
    </tr>
    <tr>
      <td>WSN-HA-plasmid</td>
      <td></td>
      <td></td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
      <td>WSN-HA-plasmid-</td>
    </tr>
  </tbody>
</table>


### Run `dms2_batch_bcsubamp`
*dms2_batch_bcsubamp* aligns the sequencing reads and counts the frequency of mutations observed.


```python
# file containing wildtype WSN HA sequence
refseq = './data/WSN-HA.fasta'

# define subamplicon alignment specifications
alignspecs = ' '.join(['1,285,36,37', 
                       '286,570,31,32',
                       '571,855,37,32',
                       '856,1140,31,36',
                       '1141,1425,29,33',
                       '1426,1698,40,43'])

# counts and alignments placed in this directory
countsdir = os.path.join(resultsdir, 'codoncounts')
if not os.path.isdir(countsdir):
    os.mkdir(countsdir)
    
# write sample information to a batch file for dms2_batch_bcsubamplicons
countsbatchfile = os.path.join(countsdir, 'batch.csv')
print("Here is the batch file that we write to CSV format to use as input:")
display(HTML(samples[['name', 'R1']].to_html(index=False)))
samples[['name', 'R1']].to_csv(countsbatchfile, index=False)

print('\nNow running dms2_batch_bcsubamp...')
log = !dms2_batch_bcsubamp \
        --batchfile {countsbatchfile} \
        --refseq {refseq} \
        --alignspecs {alignspecs} \
        --outdir {countsdir} \
        --summaryprefix summary \
        --R1trim 200 \
        --R2trim 170 \
        --fastqdir {fastqdir} \
        --ncpus {ncpus} \
        --use_existing {use_existing} 

print("Completed dms2_batch_bcsubamp.")
```

    Here is the batch file that we write to CSV format to use as input:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2B06-Lib.1</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
    </tr>
    <tr>
      <td>Mock-Lib.1</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
    </tr>
    <tr>
      <td>WSN-plasmid-</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/lgentles/1904...</td>
    </tr>
    <tr>
      <td>2B06-Lib.2</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
    </tr>
    <tr>
      <td>Mock-Lib.2</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
    </tr>
    <tr>
      <td>2B06-Lib.3</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
    </tr>
    <tr>
      <td>Mock-Lib.3</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
    </tr>
    <tr>
      <td>WSN-HA-plasmid-</td>
      <td>/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/190...</td>
    </tr>
  </tbody>
</table>


    
    Now running dms2_batch_bcsubamp...
    Completed dms2_batch_bcsubamp.


### Look at summary plots
The *_readstats.pdf summary plot shows the number of reads that were retained or thrown away due to low quality barcodes. The black regions indicate reads that failed the Illumina filter.


```python
countsplotprefix = os.path.join(countsdir, 'summary')
showPDF(countsplotprefix + '_readstats.pdf', width=700)
```


![png](analysis_notebook_files/analysis_notebook_8_0.png)


The *_readsperbc.pdf summary plot shows the number of times each barcode was read for each sample.

The distributions for this analysis look good except for the Lib.1 mock sample which could use greater sequencing depth. The Lib.2 and Lib.3 samples could also benefit from increased sequencing depth, but look pretty good for the most part.


```python
showPDF(countsplotprefix + '_readsperbc.pdf')
```


![png](analysis_notebook_files/analysis_notebook_10_0.png)


The *_bcstats.pdf summary plot shows the number of barcodes that were able to be aligned to the reference WSN HA sequence. It will also show the number of barcodes that were not able to be alligned either because the barcode was not read enough times or the sequence just could not be aligned to the reference sequence.

You can see here that the mock sample that was not sequenced to a great enough depth and has many barcodes thrown out due to too few reads.


```python
showPDF(countsplotprefix + '_bcstats.pdf', width=700)
```


![png](analysis_notebook_files/analysis_notebook_12_0.png)


The *_depth.pdf summary plot the depth (e.g. the number of reads for each codon) across the entire HA. This shows how well the barcodes were evenly represented between the different subamplicons.

The coverage looks good for most of these samples with the Lib.2 and Lib.3 having particularly even coverage.


```python
showPDF(countsplotprefix + '_depth.pdf')
```


![png](analysis_notebook_files/analysis_notebook_14_0.png)


The *_mutfreq.pdf summary plot shows the per codon frequency of mutations ar each site. In other words, it is the it is the cummulative number of mutations different from the reference sequence at each site.

We can clearly see increases in mutation frequency here in our selected samples.


```python
showPDF(countsplotprefix + '_mutfreq.pdf')
```


![png](analysis_notebook_files/analysis_notebook_16_0.png)


The *_cumulmutcounts.pdf plot below shows the fraction fo mutations that are found less than or equal to the number of times indicated on the x-axis.


```python
showPDF(countsplotprefix + '_cumulmutcounts.pdf')
```


![png](analysis_notebook_files/analysis_notebook_18_0.png)


The *_codonmuttypes.pdf summary plot the per-codon frequency of nonsynonymous, synonymous, and stop codon mutation across the entire gene.


```python
showPDF(countsplotprefix + '_codonmuttypes.pdf', width=700)
```


![png](analysis_notebook_files/analysis_notebook_20_0.png)


The *_codonmuttypes.csv summary plot shows the nummerical values for the data above.


```python
codonmuttypes = pandas.read_csv(countsplotprefix + '_codonmuttypes.csv').sort_index(axis=1)
display(HTML(codonmuttypes.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>nonsynonymous</th>
      <th>stop</th>
      <th>synonymous</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2B06-Lib.1</td>
      <td>0.001091</td>
      <td>0.000010</td>
      <td>0.000209</td>
    </tr>
    <tr>
      <td>Mock-Lib.1</td>
      <td>0.000741</td>
      <td>0.000027</td>
      <td>0.000223</td>
    </tr>
    <tr>
      <td>WSN-plasmid-</td>
      <td>0.000185</td>
      <td>0.000004</td>
      <td>0.000037</td>
    </tr>
    <tr>
      <td>2B06-Lib.2</td>
      <td>0.001183</td>
      <td>0.000011</td>
      <td>0.000248</td>
    </tr>
    <tr>
      <td>Mock-Lib.2</td>
      <td>0.000794</td>
      <td>0.000017</td>
      <td>0.000220</td>
    </tr>
    <tr>
      <td>2B06-Lib.3</td>
      <td>0.000819</td>
      <td>0.000011</td>
      <td>0.000181</td>
    </tr>
    <tr>
      <td>Mock-Lib.3</td>
      <td>0.000647</td>
      <td>0.000013</td>
      <td>0.000194</td>
    </tr>
    <tr>
      <td>WSN-HA-plasmid-</td>
      <td>0.000181</td>
      <td>0.000005</td>
      <td>0.000039</td>
    </tr>
  </tbody>
</table>


The *_codonntchanges.pdf summary plot shows the frequency of the number of nucleotide mutations that make up a given codon mutation (e.g., ATG to AAG changes 1 nucleotide, ATG to AAC changes 2 nucleotides, and ATG to CAC changes 3 nucleotides).


```python
showPDF(countsplotprefix + '_codonntchanges.pdf', width=700)
```


![png](analysis_notebook_files/analysis_notebook_24_0.png)



The *_singlentchanges.pdf plot below shows the frequency of each type of nucleotide change among only codon mutations with one nucleotide change. This plot is mostly useful to check if there is a large bias in which mutations appear. In particular, if you are getting oxidative damage (which causes G to T mutations) during the library preparation process, you will see a large excess of C to A or G to T mutations (or both). For instance, in the case of influenza, when we get bad oxidative damage, then we see an excess of C to A mutations in the final cDNA since the damage is occurring to a ssRNA genome. If you are sequencing something without polarity, you might see both types of mutations

It seems like I am seeing C to A oxidative damage in my lib 1 mock sample. This is not entirely suprising since I used older amplicons to generate this sample.


```python
showPDF(countsplotprefix + '_singlentchanges.pdf', width=700)
```


![png](analysis_notebook_files/analysis_notebook_26_0.png)


## Renumber codon counts to H3 numbering
We will now renumber the codon counts from sequential (1, 2, ... starting from the initial Met) to H3 numbering. The signal peptide is numbered with negative values, and the HA1 and HA2 subunits are numbered separately. This numbering scheme is based on an alignment to a PDB structure for an H3 HA [4HMG](https://www.rcsb.org/structure/4HMG). We will use:

1. The [./data/H1toH3_renumber.csv](./data/H1toH3_renumber.csv) file to convert the sequential numbers into the H3 numbering scheme
2. The [dms_tools2.utils.renumberSites](https://jbloomlab.github.io/dms_tools2/dms_tools2.utils.html#dms_tools2.utils.renumberSites) function to renumber the codon counts files

The renumbered files will be created in a new renumberedcounts directory and will possess the same names as the original codon counts files created from [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html) above.


```python
renumberfile = './data/H1toH3_renumber.csv'

# counts files to renumber
countsfiles = glob.glob('{0}/*codoncounts.csv'.format(countsdir))

dms_tools2.utils.renumberSites(renumberfile, countsfiles, missing='drop', 
        outdir=renumberedcountsdir)
```

## Compute differential selection
We first create a batch file to use with [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html). Note that we make the group arguments such that the samples taken from the same day with the same virus dose are compared with appropriate no antibody controls using the counts files generated above with [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html).

By grouping the replicates for the same antibody concentration and virus dose together in the batch file, we tell [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html) to analyze these together and take their mean and median.


```python
diffseldir = os.path.join(resultsdir, 'diffsel')
if not os.path.isdir(diffseldir):
    os.mkdir(diffseldir)

diffselbatchfile = os.path.join(diffseldir, 'batch.csv')
    
# create batch file for dms2_batch_diffsel
diffselbatch = pandas.DataFrame.from_records([
              # WSN library 1 selected with 50ug/mL of D045_051310_2B06
              ('2B06', '2B06-Lib.1', '2B06-Lib.1', 'Mock-Lib.1', 'WSN-plasmid-'),
              ('2B06', '2B06-Lib.2', '2B06-Lib.2', 'Mock-Lib.2', 'WSN-HA-plasmid'),
              ('2B06', '2B06-Lib.3', '2B06-Lib.3', 'Mock-Lib.3', 'WSN-HA-plasmid'),
             ],
             columns=['group', 'name', 'sel', 'mock', 'err']
             )

print("Here is the batch file that we write to CSV format to use as input:")
display(HTML(diffselbatch.to_html(index=False)))

diffselbatch.to_csv(diffselbatchfile, index=False)
```

    Here is the batch file that we write to CSV format to use as input:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>err</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2B06</td>
      <td>2B06-Lib.1</td>
      <td>2B06-Lib.1</td>
      <td>Mock-Lib.1</td>
      <td>WSN-plasmid-</td>
    </tr>
    <tr>
      <td>2B06</td>
      <td>2B06-Lib.2</td>
      <td>2B06-Lib.2</td>
      <td>Mock-Lib.2</td>
      <td>WSN-HA-plasmid</td>
    </tr>
    <tr>
      <td>2B06</td>
      <td>2B06-Lib.3</td>
      <td>2B06-Lib.3</td>
      <td>Mock-Lib.3</td>
      <td>WSN-HA-plasmid</td>
    </tr>
  </tbody>
</table>


Now we simply run [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html), getting the input counts from the directroy (--indir) where we placed the counts when we ran [dms2_batch_bcsubamp](https://jbloomlab.github.io/dms_tools2/dms2_batch_bcsubamp.html).


```python
log = !dms2_batch_diffsel \
        --summaryprefix summary \
        --batchfile {diffselbatchfile} \
        --outdir {diffseldir} \
        --indir {renumberedcountsdir} \
        --use_existing {use_existing}
```

Running this command creates a large number of output files giving the mutation and site [differential selection](https://jbloomlab.github.io/dms_tools2/diffsel.html) values in the formats of the mutdiffsel.csv and sitediffsel.csv files created by [dms2_diffsel](https://jbloomlab.github.io/dms_tools2/diffsel.html).

Specifically, there is a file for each individual sample, with a file name that gives the group and name specified for this sample in `--batchfile`.


```python
!ls {diffseldir}/
```

    2B06-2B06-50ug.log
    2B06-2B06-50ug_mutdiffsel.csv
    2B06-2B06-50ug_sitediffsel.csv
    2B06-2B06-Lib.1.log
    2B06-2B06-Lib.1_mutdiffsel.csv
    2B06-2B06-Lib.1_sitediffsel.csv
    2B06-2B06-Lib.2.log
    2B06-2B06-Lib.2_mutdiffsel.csv
    2B06-2B06-Lib.2_sitediffsel.csv
    2B06-2B06-Lib.3.log
    2B06-2B06-Lib.3_mutdiffsel.csv
    2B06-2B06-Lib.3_sitediffsel.csv
    2B06_diffsel.pdf
    2B06.log
    batch.csv
    summary_2B06-absolutesitediffselcorr.pdf
    summary_2B06-maxmutdiffselcorr.pdf
    summary_2B06-meanmutdiffsel.csv
    summary_2B06-meansitediffsel.csv
    summary_2B06-medianmutdiffsel.csv
    summary_2B06-mediansitediffsel.csv
    summary_2B06-mutdiffselcorr.pdf
    summary_2B06-positivesitediffselcorr.pdf
    summary.log
    summary_meanmaxdiffsel.pdf
    summary_meanminmaxdiffsel.pdf
    summary_meanpositivediffsel.pdf
    summary_meantotaldiffsel.pdf
    summary_medianmaxdiffsel.pdf
    summary_medianminmaxdiffsel.pdf
    summary_medianpositivediffsel.pdf
    summary_mediantotaldiffsel.pdf


Now we are going to look at the correlation of samples within a given group.

Running [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html) creates correlation plots for the mutation differential selection, positive site differential selection, absolute site differential selection, and maximum mutation differential selection for a site. These files have names like `summary_2B06-Lib.1-mutdiffselcorr.pdf`. Below we show the plots for `mutdiffsel` and `positivesiteddiffsel` (plots are also made for `absolutesitediffsel` and `maxsitediffsel`, but are not shown below as they are less informative for this experiment).



```python
diffselprefix = os.path.join(diffseldir, 'summary_')
groups = diffselbatch['group'].unique()

for seltype in ['mutdiffsel', 'positivesitediffsel']:
    print("\n{0} correlations:".format(seltype))
    plots = []
    for g in groups:
        plot = diffselprefix + g + '-' + seltype + 'corr.pdf'
        plots.append(plot)
    showPDF(plots) 
```

    
    mutdiffsel correlations:



![png](analysis_notebook_files/analysis_notebook_36_1.png)


    
    positivesitediffsel correlations:



![png](analysis_notebook_files/analysis_notebook_36_3.png)


Probably the most informative plot is simply the mean positive site differential selection. This plot show the total positive selection for all mutations combined at a give site as shown below.

This one of the most informative plots here showing strong possitive differential selection at few sites in the stem region of HA.


```python
showPDF(diffselprefix + 'meanpositivediffsel.pdf', width=800)
```


![png](analysis_notebook_files/analysis_notebook_38_0.png)


## Logo plots of differential selection
The plots above summarize the site or maximum mutation differential selection using line plots. But the most comprehensive way to show this selection is in the form of logo plots that can be created with [dms2_logoplot](https://jbloomlab.github.io/dms_tools2/dms2_logoplot.html).

We make those logog plots using the median mutation differential selection values returned by [dms2_batch_diffsel](https://jbloomlab.github.io/dms_tools2/dms2_batch_diffsel.html). The reason that we plot the median rather than the mean is that we have noticed that is is often cleaner when there are >2 replicates.

We also add underlays of the wildtype sequence. To add the wildtype sequence, we can simply specify any file with the columns `wildtype` and `site --` any of our mutdiffsel files suffice. note that we add this using the -overlay1 option, and then convert it to an underlay with `--underlay yes` since underlays look better than overlays here.

Note that this logoplot uses the converted H3 numbering asigned above.


```python
groups = diffselbatch['group'].unique()
for group in groups:
    mutdiffsel = diffselprefix + group + '-medianmutdiffsel.csv'

    logoplot = os.path.join(diffseldir, '{0}_diffsel.pdf'.format(group))
    print("\nCreating logo plot for {0} from {1}".format(group, mutdiffsel))
    
    log = !dms2_logoplot \
            --diffsel {mutdiffsel} \
            --name {group} \
            --outdir {diffseldir} \
            --restrictdiffsel positive \
            --sepline no \
            --nperline 113 \
            --overlay1 {mutdiffsel} wildtype wildtype \
            --underlay yes \
            --use_existing {use_existing}
    showPDF(logoplot)
```

    
    Creating logo plot for 2B06 from ./results/diffsel/summary_2B06-medianmutdiffsel.csv



![png](analysis_notebook_files/analysis_notebook_40_1.png)


This analysis, similarly to the positive differential selection line plot above, shows strong selection at a few sites in the stem region of HA. Sites 377, 380, 383, 384, 387, and 395 (in sequential numbering) show varing degrees of selection with site 384 appearing to have the strongest effect. Interestingly, it appears that many different mutations at this site confer escape indicating most mutations away from the wildtype sequence allow for escape. Further Dunand et al. 2015 identified escape mutations in H7 at sites V318I and I384N, which are located in the HA stalk domain, and G195E, located in the head domain. Considering the differences between H1 and H7 there are striking similarities in the selected sites. We do not identify any head escape mutations in our selection experiments here, but this could be due to strain specific differences.


```python

```
