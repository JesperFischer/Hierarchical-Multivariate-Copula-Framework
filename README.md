# Hierarchical-Multivariate-Copula-Framework

## Table of Contents
1. [Introduction](#introduction)
2. [Directory Structure](#directory-structure)
3. [Usage](#usage)


## Introduction

Here we provide all code used to generate the simulations, plots and results of the manuscript "A Hierarchical Multivariate Copula-based Framework for Cognitive Modeling".
For anonymous purposes, we have not embedded an online storage for the various models and simulation results.

## Directory Structure

The repository is structured in the following way:

```         
Thermal Pain Learning/
├── README.md             # overview of the project.
│
├── Plots/                # All code for generating figures.
│   └── ... 
│
├── Simulations/          # Directory for all simulations.
│   ├── Psychophysics/                # Directory for all simulations for the psyschophysical paradigm.
│   └──  Learning/                    # Directory for all simulations for the learning paradigm.
│
├── Real data/            # Directory containing all code for both real data set analyzed.
│   ├── Psychophysics/                # Directory for all code to fit to the the psyschophysical paradigm.
│   └──  Learning/                    # Directory for all code to fit to the the learning paradigm.
│
│
└── equations.Rmd             # Markdown for the supplementary note on equations of the models fitted.


```


## Usage

One can use the current repository to check the code for the manuscript or rerun the whole analysis pipeline. 
Upon Acceptance we will publicly share our results from the simulations and the models from the real data. 
Because rerunning the whole repository, even on a cluster, will take days.
