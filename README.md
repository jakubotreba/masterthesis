CNV Inference from scRNA-seq Workflow

Overview

This repository contains a Snakemake workflow for benchmarking four different tools for CNV inference from single-cell RNA sequencing (scRNA-seq) data. The tools included are:

SCEVAN

CopyKAT

XClone

InferCNV

The workflow processes scRNA-seq fastq files, aligns them using Cell Ranger, and generates CNV calls using the aforementioned tools.

Datasets

The workflow is tested on two datasets:

TNBC1

SNU601

