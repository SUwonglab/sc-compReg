# scCompReg #
[![license](https://img.shields.io/github/license/DAVFoundation/captain-n3m0.svg?style=flat-square)](https://github.com/DAVFoundation/captain-n3m0/blob/master/LICENSE)

scCompReg (**S**ingle-**C**ell **C**omparative **R**egulatory analysis) is an R package that provides coupled clustering and joint embedding of scRNA-seq and scATAC-seq on one sample, and performs comparative gene regulatory analysis between two conditions.

Please check the man page via `?function` (for example, `?sc_compreg`) for a detailed description of the types of inputs and outputs.

## System Requirements ##
* R (>= 3.6.0)
* Bedtools (Linux)
* Homer (Linux)

## Change Log ##
### v1.0.0 ###
* scCompReg first release.

## Intallation
Use the following command to install scCompReg R package from source code:
```R
require(devtools)
devtools::install_github("SUwonglab/sc-compReg")
```

## Full Workflow ##
The entire scCompReg workflow consists of three steps. 

The necessary inputs to `sc_compreg` are 
* consistent cluster assignments in scRNA-seq and scATAC-seq (can be obtained from coupled nonnegative matrix factorization or obtained elsewhere)
* gene expression matrices of samples 1 and 2
* chromatin accessibility matrices of samples 1 and 2
* symbol names of samples 1 and 2
* text file indicating the peak names of sample 1
* text file indicating the peak names of sample 2

In addition, download the `prior_data` directory from github via `git clone git@github.com:SUwonglab/sc-compReg.git`.


Optional: obtaining cluster assignments from coupled nonnegative matrix factorization.
1. Preproces data for `cnmf`:
    * Obtain `peak.bed` file
    * In `sc-compReg/preprocess_data/`, run the following script:
        ```bash
            bash cnmf_process_data.sh path/to/peak.bed genome_version path/to/prior_data
        ```
        where `genome_version` is one of {`hg19`, `hg38`, `mm9`, `mm10`}.
2. Run `cnmf` to get the cluster labels for sample 1 and sample 2. The cluster labels should be passed to `sc_compreg` as `O1.idx`, `E1.idx`, `O2.idx`, and `E2.idx`. For an example on how to run `cnmf`, please refer to `cnmf_example.R`



3. Obtain the rest of the necessary inputs to `sc_compreg` by editing and running in bash `sc_compreg_process_data.sh`.
4. Follow the tutorial on the `sc_compreg` function. 


## Usage ##
scCompReg provides access to the following functions:
Command       | Description
------------- | -------------
sc_compreg    | Performs single-cell comparative regulatory analysis based on scRNA-seq and scATAC-seq data from two different conditions.
mfbs_load     | Efficiently loads the `motif_target` file and returns an R `list` of the loaded objects.

## Example ##
For a full example of using the scCompReg method, please refer to `example.R`. The necessary data have been uploaded to the `data` folder in this repository.

To download the data, run the following line:
```bash
git clone git@github.com:SUwonglab/sc-compReg.git
```
The downloaded data directory will be in `sc-compReg/data/`. Simply set 
```R
path = "./sc-compReg/data/"
```
To run scCompReg, run the following lines:
```R
sample1 = readRDS(paste(path, 'sample1.rds', sep = ''))
sample2 = readRDS(paste(path, 'sample2.rds', sep = ''))

peak.name.intersect.dir = paste(path, 'PeakName_intersect.txt', sep='')
motif.target.dir = paste(path, 'MotifTarget.txt', sep='')
peak.gene.prior.dir = paste(path, 'peak_gene_prior_intersect.bed', sep='')
motif = readRDS(paste(path, 'motif.rds', sep=''))
motif.file = readRDS(paste(path, 'motif_file.rds', sep=''))

compreg.output = sc_compreg(sample1$O1,
                            sample1$E1,
                            sample1$O1.idx,
                            sample1$E1.idx,
                            sample1$symbol1,
                            sample1$peak.name1,
                            sample2$O2,
                            sample2$E2,
                            sample2$O2.idx,
                            sample2$E2.idx,
                            sample2$symbol2,
                            sample2$peak.name2,
                            motif$motif.name,
                            motif$motif.weight,
                            motif$match2,
                            motif.file,
                            peak.name.intersect.dir,
                            peak.gene.prior.dir,
                            sep.char=' ')

```
To save the obtained output, run the lines below:
```R
for (i in 1:compreg.output$n.pops) {
    write.table(compreg.output$hub.tf[[i]],
                paste(path, 'tf_', i, '.txt', sep=''),
                row.names = F,
                quote = F,
                sep = '\t')
    write.table(compreg.output$diff.net[[i]],
                paste(path, 'diff_net_', i, '.txt', sep=''),
                row.names = F,
                quote = F,
                sep = '\t')
}
```


## Citation ##
<a id="1">[1]</a> 
**Comparative regulatory analysis of single cell data reveals a novel B cell subpopulation in chronic lymphocytic leukemia**
Zhana Duren, Wenhui Sophia Lu, Joseph G. Arthur, Preyas Shah, Jingxue Xin,  Francesca Meschi, Miranda Lin Li, Corey M. Nemec, Yifeng Yin, and Wing Hung Wong