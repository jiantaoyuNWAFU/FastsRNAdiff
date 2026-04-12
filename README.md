# FastsRNAdiff

A highly efficient software for screening differentially expressed small RNA clusters.

## Table of Contents
* [Introduction](#introduction)
* [Version](#version)
* [Authors](#authors)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Usage Workflow](#usage-workflow)
* [Input Files Preparation](#input-files-preparation)
* [Command Line Options](#command-line-options)
* [Full Usage Example](#full-usage-example)
* [Output Files Explanation](#output-files-explanation)
* [Visualization Results](#visualization-results)
* [Troubleshooting](#troubleshooting)

## Introduction

FastsRNAdiff is a high-performance, user-friendly tool designed for identifying differentially expressed small RNA (sRNA) clusters between wild-type (WT) and mutant (MUT) samples. It integrates end-to-end data processing, including preprocessing (filtering invalid entries), intersection of common loci, statistical testing (Chi-square test for large samples, Fisher's exact test for small samples/zero values), q-value calculation (FDR adjustment), and multi-dimensional visualization. The tool is optimized for efficiency with parallel computing support, making it suitable for large-scale sRNA sequencing datasets.At the same time, to facilitate user adoption, we have also developed a corresponding web platform [FastsRNAdiff](http://10.105.10.204).

## Version

1.0

## Authors

Xue, Zhang and Gao etc.

Contact Email: jiantao.yu@nwafu.edu.cn

Repository URL: [https://github.com/jiaotaoyuNWAFU/FastsRNAdiff](https://github.com/jiaotaoyuNWAFU/FastsRNAdiff)

## Prerequisites

Before using FastsRNAdiff, ensure the following dependencies are installed and configured correctly.

### 1. Required Software

| Software | Version Requirement | Purpose                                | Installation Guide                                                   |
| -------- | ------------------- | -------------------------------------- | -------------------------------------------------------------------- |
| Anaconda | v25.5.1+             | Core runtime environment for the tool  | [Anaconda Official Download](https://www.anaconda.com/download/)     |
| Python   | ≥3.7                | Core runtime environment for the tool  | [Python Official Download](https://www.python.org/downloads/)        |
| cutadapt | ≥2.6                | Remove adapters or primers            | [Cutadapt Installation](https://github.com/marcelm/cutadapt) |
| samtools | ≥1.9                | Calculate Rep-total from BAM files     | [Samtools Installation](https://github.com/samtools/samtools)      |
| bowtie   | ≥1.2.3              | Align reads (for generating BAM files) | [Bowtie Installation](https://github.com/BenLangmead/bowtie) |
| bedtools | ≥2.30.0             | Unified Small RNA Cluster Coordinates  | [Bedtools Installation](http://bedtools.readthedocs.org/)
| ShortStack   | <1.0.0 (suggested)       | ShortStack:Comprehensive annotation and quantification of small RNA genes | [ShortStack](https://github.com/MikeAxtell/ShortStack) |
### 2. Required Python Libraries

FastsRNAdiff relies on the following Python libraries for data processing, statistics, and visualization. All libraries can be installed via `pip` (see [Installation](#installation) section).

| Library    | Version Requirement | Purpose                                |
| ---------- | ------------------- | -------------------------------------- |
| numpy      | ≥2.4.0             | Numerical computing (array operations) |
| scipy      | ≥1.17.0              | Statistical tests (Chi-square, Fisher) |
| pandas     | ≥3.0.0              | Data parsing and manipulation          |
| matplotlib | ≥3.10.0              | Plot generation (volcano, bar plots)   |
| seaborn    | ≥0.13.0             | Enhanced plot styling                  |
| tqdm       | ≥4.64.0             | Progress bar for long-running tasks    |
| setuptools | ≥80.0.0
| pysam      | ≥0.23.0				| Process BAM files | 

## Installation

### Command Line Runs
If you want to view the source code or run it using the command line, follow these steps to set up FastsRNAdiff on your local machine:

#### 1. Clone the Repository

First, download the tool from GitHub (replace with your actual repository URL):

```
git clone https://github.com/jiaotaoyuNWAFU/FastsRNAdiff.git
```

#### 2. Install Python Dependencies

Use the provided `requirements.txt` (see [Appendix](#installation)) to install all required libraries:

```
pip install -r requirements.txt
```

#### 3. Verify Dependencies

Confirm that all software and libraries are installed correctly:

```
# Verify Python and pip

python --version

pip --version

# Verify samtools and bowtie or ShortStack

samtools --version

bowtie --version

bedtools --version

ShortStack --version 

# Verify Python libraries

python -c "import numpy, scipy, pandas, matplotlib, seaborn, tqdm, pysam; print('All libraries installed successfully')"
```

### Web Server Runs

For the convenience of users, we provide a web-based online analysis method for downstream analysis.(see [Introduction](#introduction) section)

## Usage Workflow

The complete analysis process with FastsRNAdiff involves 5 core steps. Follow them in order to ensure valid input and accurate results:

1. **Upstream Analysis of Small RNA Sequencing Data** 

2. **Calculate Rep-total and Obtain the Counting Matrix** (use `RegionRepCalc.py`)

3. **Obtain Total Mapped Reads** (use `RegionRepCalc.py`)

4. **Downstream Analysis of Run FastsRNAdiff** (differential expression analysis)

## Input Files Preparation

FastsRNAdiff requires strictly formatted input files. Follow the steps below to prepare your data.

### 1. Step 1: Upstream Analysis of Small RNA Sequencing Data

FastsRNAdiff is well-suited for downstream analysis. For upstream analysis, users need to select the appropriate tools themselves (we recommend using ShortStack); we have also drawn inspiration from ShortStack to provide a standardized upstream analysis workflow..

#### 1.1 Remove Adapters or Primers
```
cutadapt -a <adapter_seq> [-m 18 -M 30 -q 20] -o trimmed.fastq sample1.fastq.gz(.fastq)
```
#### 1.2 Building a Reference Genome Index
```
bowtie-build genome.fa genome
```
#### 1.3 Sequence Alignment
```
bowtie [-v 1 --best --strata -k 50 -p 4] -S genome trimmed.fastq > aligned.sam
```
#### 1.4 Processing the Comparison Results
```
# 1. Convert SAM files to BAM files

samtools view -bS aligned.sam > aligned.bam

# 2. Sort

samtools sort aligned.bam -o sorted.bam

# 3. Filter out records that were not successfully mapped

samtools view -b -F 4 sorted.bam -o filtered.bam

# 4. Build index

samtools index filtered.bam
```
#### 1.5 Unify sRNA Cluster Coordinates

Standardizing the coordinates of small RNA clusters is a key prerequisite for ensuring comparability in the analysis of differential expression. If users are identifying sRNA clusters from scratch, you must merge the alignment results from multiple samples before proceeding with the identification; if users already have annotation files with standardized small RNA cluster coordinates, you may skip this step.

### 2. Step 2: Calculate Rep-total and Obtain the Counting Matrix

Rep-total (repeat-normalized total reads) is a critical metric for sRNA analysis. It normalizes read counts by accounting for multi-mapped reads (each read contributes `1/N`, where `N` is the number of mappings for that read, which can be calculated from the BAM file).

You should provide the BAM file, and use the provided `tools/RegionRepCalc.py` script to compute it.

#### 2.1 Prerequisites for the Script

* **BAM File**: Aligned reads file (The multi-mapping information must be retained.).

* **BAM Index File**: A `.bam.bai` file (generated via `samtools index <BAM>` if missing).

* **sRNA cluster File**: A tab-delimited text file with sRNA cluster coordinates. Example format:

| #Locus   |    Name[optional]     |   DicerCall[optional]        |
| -------- | ----------- | ----------------- |
| Chr1:1000-2000	|	cluster_1  | 23 |
| Chr2:3000-4000	|	cluster_2  | N |
| Chr2:10000-14000	|	cluster_3  | 21 |
| Chr3:2231019-2236105 | cluster_4 | 24 |

Among them, `#Locus` is a mandatory item, `Name` and `DicerCall` are optional. For `DicerCall`, we recommend that users adopt it as it will affect the accuracy of downstream analysis results. If users do not provide it, we will default `DicerCall` to **24nt** in the future.

#### 2.2 Run the Script

```
python RegionRepCalc.py <BAM>  <sRNAcluster> <output.txt>
```

#### 2.3 Output of the Script

A file named `<output.txt>` with the following columns:

| Column            | Description                                             |
| ----------------- | ------------------------------------------------------- |
| #Locus           | sRNA cluster coordinates (format: `ChrX:Start-Stop`)     |
| Reads      | Total number of reads aligned to the cluster            |
| Rep-total          | Repeat-normalized total reads (the value we need)    |
| UniqueReads | Count only unique mappings |
| DicerCall |  The predominant RNA sequence lengths in small RNA clusters(Valid value : **20,21,22,23,24,N,NA**) |

### 3. Step 3: Obtain Total Mapped Reads

Total mapping reads reflects the total number of reads in the sample that can be successfully aligned to the reference genome. This value can be used as a basis for subsequent standardization to correct for differences in the total number of aligned reads between different samples. When calculating this indicator, the unique mapping read segment will be fully counted, and each multiple mapping read segment will only be counted once regardless of its genome alignment frequency. This value can be automatically calculated through `RegionRepCalc.py`.

#### Key Note:

* Record Total Mapped Reads for each sample (WT and MUT).

* The order of Mapped Reads must match the order of input files (e.g., `--wt-mapped 3805767 7895813` corresponds to `--wt-files wt1.txt wt2.txt`, see Step 4 for details).

### 4. Step 4: Downstream Analysis of Run FastsRNAdiff

#### 4.1 Input Module

Both the command line and the web server must pass in the following content:

* List of wild-type sample count matrices, separated by spaces.
* Total mapped reads corresponding to the wild-type sample count matrix list in sequence.
* Mutant name(can use mutA ,mutB, etc.)
* List of mutant-type sample count matrices, separated by spaces.
* Total mapped reads corresponding to the mutant-type sample count matrix list in sequence.
* Significance level (Default: 0.05)
* Mode: Strict or non strict(Default: non strict)

Due to the software's support for simultaneously passing in different mutant types, it is necessary to distinguish them by their names.

#### 4.2 Data Preprocessing

If the input file provided by the user contains DicerCall indicators, perform pre filtering: remove small RNA clusters in each sample with DicerCall values outside the range of 20nt-24nt and Rep total greater than or equal to 1. Extract a set of common small RNA clusters with identical genome coordinates from all samples, and eliminate non common clusters from each sample.

#### 4.3 Statistical Test

Perform pairwise cross tests on samples from wild-type and mutant populations using statistical methods such as chi square test and Fisher's test. The Rep-total values of wild-type and mutant are denoted as *a* and *c*, respectively, and the total mapped reads values of wild-type and mutant are denoted as *b* and *d*, respectively. We can construct a *2 × 2* test matrix. Calculate significance level *p* and *fold change(FC)*.

| Sample group | Rep-total | Total Mapped Reads - Rep-total | 
| ------------ | ---------- | ------------------|
| wild         |  a         | b-a   |
| mutant       |  c         | d-c   |

$\log_{2}{FC} = \log_{2}{\dfrac{c/d}{a/b}}$

This software supports differential expression analysis under both strict and non strict conditions. Under non strict conditions, pre screen the collection of small RNA clusters with *p<0.05* and then perform calibration; Under strict conditions, all small RNA clusters are directly corrected. The FDR corrected value is *q-value (or adjP value)*.

#### 4.4 Differential Expression

Mark sRNA clusters with q values less than *0.05* in *m × n* independent tests as significant sRNA clusters. When $\log_{2}{FC} > 1$ in m × n test results, it is identified as a double upregulated cluster; When $\log_{2}{FC} < -1$, it is marked as a doubled down cluster; Otherwise, when $\log_{2}{FC} > 0$, it is identified as a doubled upregulated cluster; When $\log_{2}{FC} < 0$, it is marked as a doubled down cluster.

#### 4.5 Visualization

Provided volcano Plot, Bar-shaped Distribution and Sliding Windows.

## Command Line Options

FastsRNAdiff uses command-line arguments to specify input, parameters, and output. Arguments are divided into **mandatory** (required) and **optional** (customizable).

### Mandatory Parameters

These parameters must be provided; the tool will exit with an error if missing.

| Option         | Syntax                                | Description                                                                                |
| -------------- | ------------------------------------- | ------------------------------------------------------------------------------------------ |
| `--wt-files`   | `--wt-files FILE1 FILE2 ...`          | List of input files for wild-type (WT) samples (space-separated).                          |
| `--wt-mapped`  | `--wt-mapped NUM1 NUM2 ...`           | List of Mapped Reads counts for WT samples (order matches `--wt-files`).                   |
| `--mut-type`   | `--mut-type TYPE1 [--mut-type TYPE2]` | Name of mutant (MUT) type (e.g., `mutA`, `dcl2`). Use multiple times for multiple mutants. |
| `--mut-files`  | `--mut-files FILE1 FILE2 ...`         | List of input files for a MUT type (one list per `--mut-type`).                            |
| `--mut-mapped` | `--mut-mapped NUM1 NUM2 ...`          | List of Mapped Reads counts for a MUT type (order matches `--mut-files`).                  |

### Optional Parameters

These parameters override default settings (use if customization is needed).

| Option                 | Default     | Syntax                      | Description                                                            |
| ---------------------- | ----------- | --------------------------- | ---------------------------------------------------------------------- |
| `--output-dir`         | `OutputRst` | `--output-dir PATH`         | Path to store all analysis results (created automatically if missing). |
| `--significance-level` | `0.05`      | `--significance-level 0.05` | Q-value threshold for defining "significant" differential expression. |
| `--strict`          | `False`        | `--strict` | Whether to use strict mode.                           |
| `-h`/`--help`          | -           | `python FastsRNAdiff.py -h` | Show help message with all options and exit.                           |

## Full Usage Example

### Example 1: Single Mutant Type (WT vs. mutA)

```
python FastsRNAdiff.py \

	--wt-files wt_sample1.txt wt_sample2.txt \

	--wt-mapped 3805767 7895813 \

	--mut-type mutA \

    --mut-files mutA_sample1.txt mutA_sample2.txt \

    --mut-mapped 5384107 6901263 \

	--significance-level 0.05 \

	[--strict] \

    --output-dir ./AnalysisRst \
```

### Example 2: Multiple Mutant Types (WT vs. mutA + WT vs. mutB)

```
python FastsRNAdiff.py \

	--wt-files wt1.txt wt2.txt wt3.txt \

	--wt-mapped 1000000 1500000 2000000 \

	--mut-type mutA --mut-files mutA1.txt mutA2.txt --mut-mapped 900000 1100000 \

	--mut-type mutB --mut-files mutB1.txt mutB2.txt --mut-mapped 800000 950000 \

	--significance-level 0.05 \

	[--strict] \

	--output-dir ./AnalysisRst \
```

## Output Files Explanation

FastsRNAdiff organizes results into a structured directory tree for easy navigation. Below is the full structure and key file descriptions.

### Output Directory Tree

```
AnalysisRst/  # Root directory (customizable via --output-dir)

├── filterDicer/

│   ├── filterDicer_wt1_1.txt

│   └── ...

│   ├── filterDicer_mutA1_2.txt

│   └── ...

├── intersection.txt  

├── readCounts/

│   ├── readCounts_wt1_1.txt

│   └── ...

│   ├── readCounts_mutA1_2.txt

│   └── ...

├── StaticsRst/

│   ├── check_wt1_mutA1_1.txt

│   └── ...

│   ├── StaticsRst_wt_1_mutA_1.txt

│   ├── StaticsRst_wt_1_mutA_1_p_0.05.txt

│   └── ...

├── adjPvalue/

│   ├── adjPvalue_wt_1_mutA_1_p_0.05.txt

│   ├── oneColumn_parsingResult_wt_1_mutA_1_p_0.05.txt

│   └── ...

├── Separate/

│   ├── mean/

│   │   └── meanR_wt_vs_mutA.txt

│   ├── sig/

│   │   ├── sig_wt_vs_mutA.txt

│   │   └── nosig_wt_vs_mutA.txt

│   ├── up/
│   │   ├── up_meanR_wt_vs_mutA.txt

│   │   └── up2fold_meanR_wt_vs_mutA.txt

│   └── down/

│       ├── down_meanR_wt_vs_mutA.txt

│       └── down2fold_meanR_wt_vs_mutA.txt

└── Visualization/

   ├── mutA/

       ├── volcanoPlot_wt_vs_mutA.png

		├── Barshaped_wt_vs_mutA.png

       └── SlidingWindow_mutA.png
```

### Key Output Files Explained

| File Path                                     | Purpose                                                                                       |
| --------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `adjPvalue/adjPvalue_*_0.05.txt`                    | Significant differentially expressed loci (q-value < 0.05; primary result for downstream analysis). |
| `Separate/mean/meanR_wt_mutX.txt`             | Averaged metrics across biological replicates (reduces noise from single samples).            |
| `Separate/up/up2fold_meanR_*.txt`             | 2-fold up-regulated loci (stringent cutoff; high confidence results).                         |
| `Visualization/mutA/volcanoPlot_*.png` | Visually present the overall distribution characteristics of significantly differentially expressed sRNA clusters.  |
| `Visualization/mutA/Barshaped_*.png` | Visually present the distribution and expression level of significantly differentially expressed sRNA clusters on chromosomes.  |
| `Visualization/mutA/SlidingWindow_*.png` | Characterize the relative density of sRNA clusters within the current sliding window.                   |

### 1. Volcano Plot

* **Purpose**: Show the relationship between fold change (x-axis) and statistical significance (-log10(q-value); y-axis).

* **Guides**:

  * Horizontal dashed line: q-value = 0.05 (-log10(0.05) ≈ 1.3)

  * Vertical dashed lines: $\log_{2}{FC} = ±1$ (2-fold change cutoff)

### 2. Bar-shaped Distribution

* **Purpose**: Map differentially expressed loci to their chromosomal positions and show regulation direction (up/down).

* **X-axis**: Chromosomal position (in Mb, megabases).

* **Y-axis**: log2FC (positive = up-regulated, negative = down-regulated).

### 3. Sliding Window

* **Purpose**: Reduce clutter in large datasets by counting significant loci in sliding windows (avoids overplotting).

* **Parameters**: Window size = 1 Mb, step size = 0.1 Mb (fixed).

* **X-axis**: Chromosomal position (Mb).

* **Y-axis**: Normalized density of significant loci (0 to 1; higher = more loci in the window).

## Troubleshooting

### Debugging Tips

1. **Check Input Formats**: Verify that all input files have the correct header (e.g., `#Locus` instead of `Locus`; tab-delimited, not comma-delimited).

2. **Inspect Intermediate Files**:

* `intersection.txt`: If empty, no common loci exist across samples (check if input files are from the same genome build).

* `StaticsRst/StaticsResult_*.txt`: If missing, the statistical test failed (check for non-numeric values in `readCounts/` files).

1. **Contact Support**: For persistent issues, send the following to jiantao.yu@nwafu.edu.cn:

* Full command line used.

* Error message (copy-paste the entire traceback).

* Header and 2-3 data lines of your input files.




