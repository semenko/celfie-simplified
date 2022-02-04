
# CelFiE — Refactored
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/semenko/celfie-simplified/blob/master/notebooks/celfie-simplified-demo.ipynb)


**By Nick Semenkovich \<semenko@alum.mit.edu\>**

<br>

An opinionated rework of the CelFiE project (original repository [here](https://github.com/christacaggiano/celfie)). Note that this project simply reworks CelFiE to make it more Pythonic / reproducible and accept standard (.bed) inputs.

The goal of this code is to predict fractional tissue abundance from a mixed population of cells, using cell free methylation data. You likely need to build your own reference matrix, detailed below, though the original is available for reference.

## Usage


This repo can be run entirely in Google Colab: [celfie-simplified-demo.ipynb](https://colab.research.google.com/github/semenko/celfie-simplified/blob/master/notebooks/celfie-simplified-demo.ipynb)


Locally, you can run:
```
git clone https://github.com/semenko/celfie-simplified
cd celfie-simplified

pip3 install -r requirements.txt

python3 celfie-simplified.py --input data/sample-neutrophil.bed --reference_tims data/caggiano_TIM_matrix.bed --unknowns 0 --output sample-output/sample-neutrophil
```

## Input Formats

CelFiE expect two input files: your sample's methylation data (as a .bed) and the tissue informative marker (TIM) matrix.

Your sample's data .bed should have columns 4 and 5 set to the **# of methylated reads** and # of total reads:

```
# chr start end   Hepatocyte_meth Hepatocyte_depth
chr1	10    11	44	63
chr1	50    51	71	133
chr1	60    61	89	115
```

**Note:** Your sample .bed does not need a header. If you do not provide one, samples will be named "sample1, sample2 …". If you provide one, it **must** start with #, and your sample names must be formatted as "Tissue_name_meth" and "Tissue_name_depth". You can include more than one tissue (e.g. columns 5 and 6 can be tissue2_meth and tissue2_depth).

**Note**: Many analyses generate % methylation values — you can convert from percent to absolute counts awk:
```
XXXX
```

## TIM Matrix Format

CelFiE expects tissue informative markers (TIMs) in a .bed file with a header the following format:

```
# chrom start  end      tissue1_meth      tissue1_depth     tissue2_meth      tissue2_depth
chr1	10	11	25	29    53	105
chr1	50	51	85	99    72	285
chr1	60	61	92	117   12	33
```

**Note:** Here, the .bed header is not optional, as CelFiE needs to have valid tissue names for assignments and subsequent plotting.

As a reference, you can use the TIM matrix from the [original CelFiE paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5270101/), which is in `data/caggiano_TIM_matrix.bed`, but this matrix has some caveats (see below).


### Outputs

Output formatting is unchanged from the original CelFiE code. CelFiE outputs tissue estimates for each sample in your input — - i.e. the proportion of each tissue in the reference making up the cfDNA sample. See `celfie_demo/sample_output/1_tissue_proportions.txt` for an example of this output.

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

## Custom TIM Matrix Generation

The original TIM matrix is in `data/ was trained on a combination of human hg38 data from ENCODE and Blueprint, as described in the original paper. It encompasses 19 human tissues:

However, this matrix is ***




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

Please note all data was converted to hg38 and all CpGs are reported as (Chrom, start, end), where the end position indicates the C in the CpG dinucleotide.  

#### Selecting TIMs

Code to find TIMs is located at `TIMs/tim.py`. This code takes a reference bedfile of all the tissues you would like to calculate TIMs for as input. See `TIMs/sample_input.txt.`

The TIM code can be run as:

```bash
python tim.py <input file> <output file> <num of tim/tissue> <num of tissues> <depth filter> <nan filter>
```

The number of TIMs per tissue can be adjusted, but note that as the number of TIMs approaches the number of CpGs, the less informative that TIM will be for that tissue.

The **depth filter** only will consider CpGs that have a median depth across all tissues greater than a user specified value. This is to ensure that low-coverage CpGs do not get selected as TIMs. The **NaN filter** will only consider CpGs that have less than a user specified number of missing values. This is to ensure a TIM isn't selected for a tissue because it is one of the few tissues with data at that location. The **number of tims/tissue** can vary. We find that 100 is a good number, and note that as the number of TIMs increase, the lower quality the TIMs will be, since we are selecting the top most informative CpGs/tissue (in other words, the top 100 most informative CpGs for pancreas will by definition, be "better" than the top 500).


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
