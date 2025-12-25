# Model-based predictions of protective HIV PrEP levels

>This repository contains the code used to analyze the prophylactic efficacy of TDF/FTC in cisgender women and men who have sex with men (MSM). It supports two analytical strategies: the top-down and bottom-up approaches.
>
> For cisgender women: details described in the associated paper (DOI: 10.21203/rs.3.rs-2772765/v1)
> 
> For MSM: the top-down and bottom-up analysis are applied on MSM. A new top-down approach is described in the submitted manuscript.
>
>The code provided enables full reproduction of the top-down and bottom-up analyses, and the data contained can be used to reproduce the figures in the corresponding paper.

[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) [![DOI](https://zenodo.org/badge/616463573.svg)](https://zenodo.org/badge/latestdoi/616463573)

## Table of Contents
-   [System requirements](#system-requirements)
      -   [Operating systems](#operating-systems)
      -   [Prerequisites](#prerequisites)
      -   [Dependencies](#dependencies)
- [PrEP MSM](#PrEP_MSM)
- [PrEP women](#PrEP_women)

## System requirements

### Operating systems
This workflow was tested on Ubuntu 20.04.5 LTS and Mac OS.

### Prerequisites
Some tools have to be installed to run the analysis. We recommend following the steps outlined below.

#### Install Conda/Miniconda

Conda will manage the dependencies of our program. Instructions can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install.


#### Create the working environment

Create a new environment from the given environment config in [`env.yml`](./env/env.yml), where the pipeline will be executed.
Go to the main folder of the repository and call:

```
conda env create -f env/env.yml
```

This step may take a few minutes.

To activate the environment type:

```
conda activate prep
```

### Dependencies

This workflow uses the following dependencies:

```
mpmath==1.3.0
numpy==2.4.0
pandas==2.3.3
scipy==1.16.3
torch==2.2.0
```
They can be installed automatically by creating the conda environment above. 

## PrEP MSM

Contains all codes to reproduce analysis in the manuscript "Pharmacological markers of HIV prevention for oral pre-exposure prophylaxis in MSM" by Sara Iannuzzi, Malin Müller, Yifan Yu, Lanxin Zhang, Craig W. Hendrix, Robert R.
Bies and Max von Kleist (https://doi.org/10.21203/rs.3.rs-7537396/v1).

## PrEP women

Contains all codes to reproduce analysis in the manuscript "Model-based predictions of protective HIV pre-exposure prophylaxis adherence levels in cisgender women" by Lanxin Zhang*, Sara Iannuzzi*, Ayyappa Chaturvedula, Elizabeth Irungu, Jessica E Haberer, Craig W. Hendrix and Max von Kleist, Nat Med. 2023 Nov;29(11):2753-2762 (preprint; https://doi.org/10.1038/s41591-023-02615-x).
