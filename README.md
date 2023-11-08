# TeloSeq

TeloSeq is a pipeline for the analysis of telomere sequencing data generate ONT on R9 and R10 pore.
A custom model was trained on the R9 data to improve the accuracy of the basecalling telomeres.

## Requirements
- conda>=23.5.0
- python>=3.6
- [NCRF](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder)
- ONT [Bionto basecaller](https://github.com/nanoporetech/bonito) == 0.6.2

## Installation

### Noise cancelling Repeat Finder (NCRF)

```
add the path to NCRF to your .bashrc file
```

### TeloSeq pipeline


```
git clone https://github.com/Priyesh000/teloseq.git
cd teloseq
conda env create -f environment.yml 
conda activate teloseq
```
