# MSI-H Colorectal Cancer Frameshift Neoantigen Prediction

A self-directed introductory bioinformatics pipeline to identify and rank candidate MHC Class I frameshift neoantigens from microsatellite-instable (MSI-H) colorectal cancer mutation data.

## Overview
- **Input Data:** Frameshift insertion/deletion mutations from TCGA-COAD MSI-H samples.
- **Sequence Translation & Filtering:** Translated wild-type and mutated coding sequences via Ensembl REST API and Biopython, extracting **584 novel frameshift 9-mers** absent from wild-type proteins.
- **MHC Class I Binding Prediction:** Scored peptide binding affinities across four HLA Class I alleles (`HLA-A*02:01`, `HLA-A*11:01`, `HLA-A*24:02`, `HLA-B*07:02`) using **MHCflurry 2.3.3** (2,016 peptide–allele combinations).
- **Top Candidates:** Identified **75 candidate binders** (`<= 500 nM`), including high-affinity hits (`< 50 nM`) in `ENTPD2`, `SOX9`, `KMT2D`, and `LARP4B`.

## Key Files
- `generate_9mers_fast.py` — Fetches canonical CDS, simulates frameshift translation, and filters novel 9-mers.
- `predict_mhc.py` — Runs MHCflurry affinity predictions across selected HLA alleles.
- `top_neoantigen_candidates.csv` — Ranked table of predicted binders (`<= 500 nM`).
- `neoantigen_summary.png` — Summary visualization of top candidate 9-mers and gene distribution.
-# msih-colorectal-neoantigen-prediction
self-directed bioinformatics pipeline, predicting MHC Class I frameshift neoantigens in MSI-H colorectal cancer
