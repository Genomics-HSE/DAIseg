## DAIseg

### Overview 

DAIseg infers introgressed archaic segments in modern genomes.

### Citation

Planche, L., Ilina, A.V., & Shchur, V.L. (2024). Highly Accurate Method for Detecting Archaic Segments in the Modern Genomes. *Lobachevskii J Math*, 45, 2910–2917. DOI: [10.1134/S1995080224602959](https://doi.org/10.1134/S1995080224602959).

### Installation

DAISeg is a Python package, so you can install it from Github to your Python environment.

We recommend using [`pixi`](https://pixi.prefix.dev/latest/installation/] to manage your research environment. With `pixi`, DAIseg can be installed with 

``` bash
pixi init
pixi add --git https://github.com/Genomics-HSE/DAIseg --pypi daiseg
```

Beyond Python dependencies, DAIseg needs `bedtools` and `bcftools` (for data preprocessing steps only).
Example `pixi` environment with all required dependencies is provided in the `examples/` directory, see "Usage" below.

### Usage 

Here is an example of how to run DAIseg on chromosome 22 with hg19 1000Genomes samples.

**0.** Install `DAIseg` in an environment of your choice. Here is how I would do it. 

First, make sure you have [`github cli`](https://github.com/cli/cli#installation) and [`pixi`](https://pixi.prefix.dev/latest/installation/) installed on your Linux machine.

Clone this repository. Then, copy the example folder out of the git repository to use it as your working directory (let's call it `example_run`).

``` bash
gh repo clone Genomics-HSE/DAIseg
cp -r DAISeg/examples/grch37 example_run
cd example_run
```

Install the `pixi` environment.

``` bash
pixi install
```

Test the environment, making sure DAIseg launches.

``` bash
pixi run daiseg
```

This command should do nothing by design. Email us if you have any problems.



**1.** Download the data and fill out the config file. Run the `download.sh` script to get all the necessary data. It will create the following directory layout:

```
$ tree data
data
├── 1000g
│   ├── 20140520.strict_mask.autosomes.bed
│   ├── 20140520.strict_mask.chr22.bed
│   ├── ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
│   └── ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
├── ancestral
│   ├── homo_sapiens_ancestor_22.bed
│   ├── homo_sapiens_ancestor_22.fa
│   └── ....
├── gaps
│   └── gap.txt
├── nd
│   ├── altai
│   │   ├── chr22_mask.bed.gz
│   │   ├── chr22_mq25_mapab100.vcf.gz
│   │   └── chr22_mq25_mapab100.vcf.gz.tbi
│   └── vindija
│       ├── chr22_mask.bed.gz
│       ├── chr22_mq25_mapab100.vcf.gz
│       └── chr22_mq25_mapab100.vcf.gz.tbi
├── reference
│   └── hs37d5.fa.gz
└── sizes
    └── hg19.chrom.sizes
```

The folder will be ~5.2GiB in size. 

**1.** Enter data file locations and inference parameters into a `json` file. The provided `example_config.json` file already matches the `data` directory created above. The inference parameters in the file list the YRI population samples in the `outgroup` field and IBS population samples in the `ingroup` field. See below for the specification of config fields. 

**2.** Extract the relevant samples from the 1000Genomes files. Run 
```bash
pixi run daiseg restrict_1kG -json example_config.json -threads 8
```

This will create the `<prefix>/preprocessed.vcf.gz` file in your output prefix directory (specified in the config JSON, see below). This is the filtered file that contains only the modern human samples listed as `outgroup` and `ingroup` in the config


**3.** Create callability masks. Run 

``` bash
pixi run daiseg callability -json example_config.json
```

This will create the callability/coverage for modern samples (`<prefix>/1kg_coverage.bed`) and for Neanderthal samples (`<prefix>/nd_1kg_coverage.bed`).

**4.** Prepare input data for the HMM. Run 

``` bash
pixi run daiseg prep -json example_config.json -threads 8
```

This will create `<prefix>/preprocessed_data.tsv` -- the file tabulating all sites that will be used by the HMM.

**5.** Run the algorithm. There are three available options:

- `daiseg run -json <config>` will do .... This command takes a single JSON file as an argument.

- `daiseg run_EM -jsons <configs>` will do .... This command can be used on a batch of config JSON files, so that multiple chromosomes can be processed at once.

- `daiseg run_EM_trans -jsons <configs>` will .... This commans also can be used with a batch of config JSONs.

Run the EM algorithm:

``` bash
pixi run daiseg run_EM -jsons example_config.json
```

The output file `<prefix>/out.hg19.chr22.em.tsv` will contain inferred archaic segments for each ingroup sample.

### Config specification

- `CHROM`: chromosome name
- `prefix`: output location, relative to the working directory; all intermediate and output files will be in this directory
- `output`: name prefix for output files containing inferred segments; will be put under the output prefix directory
- `data`: name for the HMM algorithm input data file; will be created automatically by preparation steps and put in the output prefix directory
- `window_callability`: filepaths for window vallability `bed` files; these files will be created by the preprocessing steps, so there is no need to change these
- `files`: filepaths for input data files; see the `example/` directory for expected file formats and download links for typically used files
  * `neand_files`: `.vcf` variant file and `.bed` accessibility mask file for each neanderthal genome; see example config in `examples/` for the format
  * `1000GP_files`: `bed` is the location of the 1000Genomes strict accessibility mask; `vcf_initial` is the location of the 1000Genomes variant file for the focal chromosome; `vcf` is the filename for the preprocessed `vcf`, it will be created automatically and put under the output prefix so there is usually no need to change this field
  * `ancestral.fasta`: fasta file for _Homo sapiens_ ancestor 
  * `reference.fasta`: reference _Homo sapiens_ fasta file
  * `chr_lengths`: chromosome lengths
- `gaps`: filepath for gap annotations for the genome assembly in use; typically downloaded from the UCSC database, see `download.sh` script in `examples/`
- `samples`: modern hyman sample names to use in the analysis
  * `outgroup`: list of sample names found in the 1000Genomes vcf to use as outgroup (population with no archaic admixture)
  * `ingroup`: list of sample names to use as ingroup (population with archaic admixture)
  * `neand`: list of sample names in the Neanderthal `vcf` files
- `parameters_initial`: demographic parameters of the model, with times listed in years

