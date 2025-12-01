# Malaria Host Response Transcriptomic Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.xxxx/zenodo.xxxxx.svg)](https://doi.org/10.xxxx/zenodo.xxxxx)

## 📋 Overview

Comprehensive transcriptomic analysis of host response to *Plasmodium falciparum* infection. This study identifies key immune pathways, differential gene expression patterns, and potential biomarkers associated with malaria pathogenesis.

### Key Findings:
- Identified 1,247 differentially expressed genes in malaria-infected samples
- Discovered enrichment in interferon signaling, cytokine production, and inflammatory response pathways
- Revealed novel host factors potentially involved in malaria severity

## 🎯 Objectives

1. Characterize host transcriptional response to malaria infection
2. Identify differentially expressed genes between infected and control samples
3. Perform functional enrichment analysis (GO, KEGG, Reactome)
4. Conduct Gene Set Enrichment Analysis (GSEA) using MSigDB collections
5. Identify potential biomarkers and therapeutic targets

## 📁 Project Structure

malaria-transcriptomic-analysis/
│
├── data/
│   ├── raw/                    
│   │   ├── GEO_metadata.txt    # Sample information
│   │   └── README_data.md      # Data provenance and sources
│   ├── processed/              # Processed data
│   │   ├── normalized_counts.csv
│   │   ├── deg_results.csv
│   │   └── sample_metadata.csv
│   └── external/               # Reference files
│       ├── gene_annotations.csv
│       └── pathway_databases/
│
├── code/
│   ├── 01_data-preprocessing.R
│   ├── 02_quality-control.R
│   ├── 03_differential-gene-expression-analysis.R
│   ├── 04_functional-enrichment-analysis.R
│   ├── 05_pathway_analysis.R
│   ├── 06_gsea_analysis.R
│   ├── 07_visualization.R
│   └── run_analysis.sh         # Pipeline script
│
├── results/
│   ├── figures/                # Publication-quality figures
│   │   ├── qc_plots/
│   │   ├── volcano_plots/
│   │   ├── heatmaps/
│   │   ├── enrichment_plots/
│   │   └── pathway_analysis/
│   ├── tables/                 # Analysis outputs
│   │   ├── differential_expression/
│   │   ├── go_enrichment/
│   │   ├── kegg_pathways/
│   │   └── gsea_results/
│   └── reports/                # Generated reports
│       └── exploratory_analysis.html
│
├── docs/                       # Documentation
│   ├── methods.md
│   ├── data_sources.md
│   ├── interpretation.md
│   └── references.bib
│
├── notebooks/                  # Jupyter/R Markdown notebooks
│   ├── exploratory_analysis.Rmd
│   └── final_analysis.Rmd
│
├── environment/
│   ├── environment.yml         # Conda environment
│   ├── Dockerfile              # Container specification
│   └── sessionInfo.txt         # R session info
│
├── tests/                      # Unit tests for functions
│   └── test_analysis_functions.R
│
├── LICENSE
├── README.md                   # Main repository documentation
├── CITATION.cff                # Citation file
└── .gitignore
