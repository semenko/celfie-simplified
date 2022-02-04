
# CelFiE — Refactored
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/semenko/celfie-simplified/blob/master/notebooks/celfie-simplified-demo.ipynb)


**By Nick Semenkovich \<semenko@alum.mit.edu\>**

<br>

An opinionated rework of the CelFiE project, original repository [here](https://github.com/christacaggiano/celfie).

The goal of this code is to determine fractional tisuses in a population of cells, from cell free methylation data. This adapts the CelFiE code to take standard inputs (.bed), and guts a lot of code. 

This only implements CelFiE's EM algorithm and the tissue informative marker (TIM) generation.

# Usage


This repo can be run entirely in Google Colab: [celfie-simplified-demo.ipynb](https://colab.research.google.com/github/semenko/celfie-simplified/blob/master/notebooks/celfie-simplified-demo.ipynb)

## Input Data Format

CelFiE expects the methylation data to be in the form # of methylated reads, # of total reads. For example it could look like:

```
CHR   START END METH DEPTH
chr1	10	11	44.0	63.0
chr1	50	51	71.0	133.0
chr1	60	61	89.0	115.0
```

## TIM Matrix Format

CelFiE should work, in theory, on Illumina Chip data, if you estimate the read depth of each of the sites. However, we do not officially recommend this.

The input of CelFiE is a single txt file including both the reference data and the cfDNA, with a header indicating sample names (see `celfie_demo/sample_data.txt`). Essentially the file is the reference and cfDNA sample bed files combined. This data should look something like this:


```
CHROM START END SAMPLE1_METH SAMPLE1_DEPTH CHROM START END TISSUE1_METH TISSUE1_DEPTH
chr1	10	11	44.0	63.0  chr1	10	11	25.0	29.0
chr1	50	51	71.0	133.0 chr1	50	51	85.0	99.0
chr1	60	61	89.0	115.0 chr1	60	61	92.0	117.0
```


### Output

CelFiE will output the tissue estimates for each sample in your input - i.e. the proportion of each tissue in the reference making up the cfDNA sample. See `celfie_demo/sample_output/1_tissue_proportions.txt` for an example of this output.

```
        tissue1 tissue2 .... unknown
sample1 0.05 0.08 .... 0.1
sample2 0.7 0.12 .... 0.2

```

CelFiE also outputs the methylation proportions for each of the tissues plus however many unknowns were estimated. This output will look like this:

```   
      tissue1  tissue2 ... unknown
CpG1  0.99 1.0 ... 0.3
CpG2  0.45 0.88 ... 0.1
```

Sample code for processing both of these outputs can be seen in `demo.ipynb`.

## Tissue Informative Markers

TIMs are available at `TIMs/sample_tims.txt` for individual CpG TIMs, and `TIMs/sample_tims_summed.txt` for reads summed +/-250bp around a TIM. We recommend using the `TIMs/sample_tims_summed.txt` for improved decomposition performance.

The TIMs represent markers for the following tissues:

- dendritic cells
- endothelial cells
- eosinophils
- erythroblasts
- macrophages
- monocytes
- neutrophils
- placenta
- T-cells
- adipose
- brain
- fibroblasts
- heart left ventricle
- hepatocytes
- lung
- mammary gland
- megakaryocytes
- skeletal muscle myoblasts
- small intestine

Data was retrieved from the [ENCODE](https://www.encodeproject.org/search/?type=Experiment&assay_slims=DNA+methylation&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=WGBS) and [Blueprint](http://dcc.blueprint-epigenome.eu/#/files) data portals. When available, two biological replicates per tissue were combined into one sample. The TIMs were then calculated on the combined sample.

Please note all data was converted to hg38 and all CpGs are reported as (Chrom, start, end), where the end position indicates the C in the CpG dinucleotide.  

#### Selecting TIMs

Code to find TIMs is located at `TIMs/tim.py`. This code takes a reference bedfile of all the tissues you would like to calculate TIMs for as input. See `TIMs/sample_input.txt.`

The TIM code can be run as:

```bash
python tim.py <input file> <output file> <num of tim/tissue> <num of tissues> <depth filter> <nan filter>
```

The number of TIMs per tissue can be adjusted, but note that as the number of TIMs approaches the number of CpGs, the less informative that TIM will be for that tissue.

The **depth filter** only will consider CpGs that have a median depth across all tissues greater than a user specified value. This is to ensure that low-coverage CpGs do not get selected as TIMs. The **NaN filter** will only consider CpGs that have less than a user specified number of missing values. This is to ensure a TIM isn't selected for a tissue because it is one of the few tissues with data at that location. The **number of tims/tissue** can vary. We find that 100 is a good number, and note that as the number of TIMs increase, the lower quality the TIMs will be, since we are selecting the top most informative CpGs/tissue (in other words, the top 100 most informative CpGs for pancreas will by definition, be "better" than the top 500).

For the sample data provided, we suggest:

```bash
python tim.py sample_input.txt tim.txt 100 19 15 2
```

#### Combining Reads

In our paper, we found that summing all reads +/-250bp offered improved performance when decomposing. To do this for TIMs generated as output of `tim.py`, we provide a shell script `TIMs/tim.sh` to call TIMs and sum data.

This script can be updated to change the following parameters:

```bash
input_file=sample_input.txt
output_file=sample_tims.txt
summed_file=sample_tims_summed.txt
sum_window_size=500
number_tims=100
number_tissues=19
depth_filter=15
na_filter=2
```

The pipeline can then be ran as
```bash
./tim.sh
```

## Citation

Christa Caggiano, Barbara Celona, Fleur Garton, Joel Mefford, Brian Black, Catherine Lomen-Hoerth, Andrew Dahl, Noah Zaitlen, *"Comprehensive cell type decomposition of circulating cell-free DNA with CelFiE"*, Nature Communications, May 2021 [link](https://www.nature.com/articles/s41467-021-22901-x)
