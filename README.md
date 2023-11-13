# TeloSeq

TeloSeq is a comprehensive pipeline designed for the in-depth analysis of telomere sequencing data generated using ONT on both R9 and R10 pores. To enhance basecalling accuracy for telomeres, a custom model has been trained for R9 Pores. 

The pipeline operates on basecalled fastq files and progresses through multiple stages, culminating in the generation of critical outputs such as the telomere length distribution, mapped telomere reads, Telomere Variant Repeat (TVR) signal, and various statistical metrics.

It is imperative to note that, specifically for R9 data, the utilization of the custom model is indispensable during the basecalling step. Conversely, for R10 data, default basecalling using the SUP model as an integral component of the Dorado basecaller suffices. This nuanced approach ensures optimal performance and accuracy tailored to the characteristics of each sequencing dataset.


## ADD CONTENTS TABLE HERE

TODO:
- [ ] Add contents table
- [ ] Add installation instructions
- [ ] Add usage instructions
- [ ] Add example config file
- [ ] Add example scripts
- [ ] Add example data
- [ ] Add example output


# System requirements

## Hardware requirements
TeloSeq pipeline requires only a standard computer with sufficient RAM to support the in-memory operations. For minimal performance, this should be at least 8 GB of RAM. For optimal performance, we recommend 16 GB or more of RAM. The pipeline is compatible with multi-core systems. However, the pipeline is not optimized for parallel computing.

## Software requirements
### OS requirements
This pipeline has been tested on Linux operating systems (Ubuntu 16.04+). The pipeline is not currently supported on Windows.

### Other Requirements
- conda>=23.5.0
- python>=3.6
- [NCRF](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder) == v1.01.00
- ONT [Bionto basecaller](https://github.com/nanoporetech/bonito) == 0.6.2, if using R9.4.1 data, or [Dorado basecaller](https://github.com/nanoporetech/dorado) >= 0.4.0, if using R10 data


## Installation guide

### TeloSeq pipeline

```
git clone https://github.com/Priyesh000/teloseq.git
cd teloseq
conda env create -f environment.yml 
conda activate teloseq
```

### Noise cancelling Repeat Finder (NCRF)

To install NCRF from source, utilize the provided command. Further details and instructions can be accessed [here](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder).
After the installation process is finished, it is essential to append the file path to the [config.yml](config.yml) file. Locate the  [config.yml](config.yml) file at ./config.yml.

```
git clone --branch v1.01.00 https://github.com/makovalab-psu/NoiseCancellingRepeatFinder.git  
cd NoiseCancellingRepeatFinder  
make  
```

## Getting started

### Build reference genome
To build the reference genome, use the following command. This command will generate a directory named "refgenome" within the current working directory. The reference genome is constructed using data from the Telomere-to-Telomere consortium (T2T), specifically utilizing the HG002v0.7 and T2T CHM13 reference genomes. This process results in the creation of two distinct files: one for the complete reference genome and another for the paternal reference genome. The paternal reference genome is particularly valuable for telomere phasing. The full reference genome, on the other hand, is employed for read mapping and the assignment of telomeres to chromosomes. Each reference genome is accompanied by a TSV file that contains information about the start and end positions of the telomeres within the reference genome. Additionally, the chromosomes are further divided into three sections: P arms, Q arms, and the M (middle) region.


After the references are built, you should add the file paths to the [refgenome.csv](config/refgenome.csv) and provide unique names (`refgenome_id`) for each reference genome. These unique names (`refgenome_id`) will be used in the [datasets.csv](datasets.csv) file to specify which reference genome should be utilized for the analysis. You can find the refgenome.csv file at config/refgenome.csv.

```

python scripts/build_telomere_reference.py -o refgenome/ full

python scripts/build_telomere_reference.py -p -o refgenome/ hg002_paternal

```

**example refgenome.csv file**

```
refgenome_id,refgenome_path,anchor_tsv
full,refgenome/full.fasta,refgenome/full.tsv
hg002_paternal,refgenome/HG002_paternal.fasta,refgenome/HG002_paternal.tsv
```    


### Dataset configuration
Here is where we will specify the path to our base called data and the reference genome to be used for the analysis. The [datasets.csv](dataset.csv) file is located at [./datasets.csv](dataset.csv).
Please note that the `refgenome_id` must match the `refgenome_id` in the [refgenome.csv](config/refgenome.csv) file.


### Run TeloSeq pipeline
Teloseq pipeline can be run using the following command. The output data will to saved at the location given in [config.yml](config.yml) under `SAVEPATH`. The output directory will contain the final outputs of the pipeline under `SAVEPATH/reports/*/*.html` directroy. The final outputs include the telomere length distribution, mapped telomere reads, TVR signal, and various statistical metrics. The pipeline will also generate intermediate files that can be used for further analysis. The intermediate files include the telomere length distribution for each chromosome, the TVR signal for each chromosome, and the mapped telomere reads for each chromosome. The intermediate files are located in the `SAVEPATH/*` directory. The pipeline will also generate a log file that contains information about the pipeline run. The log file is located in the `.snakemake/logs/*` directory.

```
snakemake --use-conda --configfile config.yml  -pr --cores 8 
```
