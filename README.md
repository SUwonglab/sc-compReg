# scCompReg #
[![license](https://img.shields.io/github/license/DAVFoundation/captain-n3m0.svg?style=flat-square)](https://github.com/DAVFoundation/captain-n3m0/blob/master/LICENSE)

scCompReg (**S**ingle-**C**ell **C**omparative **R**egulatory analysis) is an R package that provides coupled clustering and joint embedding of scRNA-seq and scATAC-seq on one sample, and performs comparative gene regulatory analysis between two conditions.

Please check the man page via `<?function>` (for example, `<?sc_compreg>`) for a detailed description of the types of inputs and outputs.

## System Requirements ##
* R (>= 3.6.0)

## Change Log ##
### v1.0.0 ###
* scCompReg first release.

## Intallation
Use the following command to install scCompReg R package from source code:
```R
require(devtools)
devtools::install_github("SUwonglab/sc-compReg")
```

## Usage ##
scCompReg provides access to the following functions:
Command       | Description
------------- | -------------
sc_compreg    | Performs single-cell comparative regulatory analysis based on scRNA-seq and scATAC-seq data from two different conditions.
mfbs_load     | Efficiently loads the `<motif_target>` file and returns an R `<list>` of the loaded objects.

## Example ##
For an example of using the scCompReg method, please refer to `<example.R>`. The necessary data have been uploaded to the `<data>` folder in this repository.


## Citation ##
<a id="1">[1]</a> 
**Comparative regulatory analysis of single cell data reveals a novel B cell subpopulation in chronic lymphocytic leukemia**
Zhana Duren, Wenhui Sophia Lu, Joseph G. Arthur, Preyas Shah, Jingxue Xin,  Francesca Meschi, Miranda Lin Li, Corey M. Nemec, Yifeng Yin, and Wing Hung Wong