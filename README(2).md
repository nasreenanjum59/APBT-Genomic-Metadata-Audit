# APBT Genomic Metadata Audit

Reproducibility code and data for the study:

**From data poisoning to advanced persistent biological threats (APBTs): a taxonomy and large-scale metadata audit of genomic surveillance infrastructure**

This repository contains the analysis programs used to retrieve, audit, enrich, and analyse SARS-CoV-2 BioSample metadata within the Advanced Persistent Biological Threats (APBT) framework.

## Repository structure

```text
APBT-Genomic-Metadata-Audit/
├── code/
│   ├── 01_APBT_BioSample_Metadata_Retrieval.py
│   ├── 02_APBT_BioSample_Local_Audit.py
│   ├── 03_APBT_CrossRepository_Linkage_Enrichment.py
│   └── 04_APBT_CrossRepository_Final_Analysis.py
│
└── data/
    ├── README.md
    └── SARS_CoV2_BioSample_metadata_100000.csv
```

## Programs

### Program 1: BioSample metadata retrieval

`code/01_APBT_BioSample_Metadata_Retrieval.py`

Retrieves SARS-CoV-2 BioSample metadata from NCBI using the Entrez API and writes the parsed records to CSV.

### Program 2: BioSample-local audit

`code/02_APBT_BioSample_Local_Audit.py`

Runs the 23-check BioSample-local metadata audit across five metadata layers. This program is retained as a methodological baseline and does **not** represent the final cross-repository results reported in the revised manuscript.

### Program 3: Cross-repository linkage and enrichment

`code/03_APBT_CrossRepository_Linkage_Enrichment.py`

Extends selected provenance and technical checks to linked NCBI BioProject, SRA, and Nucleotide/GenBank resources. The program records whether information was present directly in BioSample, recovered from a linked resource, or not found in the queried resources.

The status `not_found_in_queried_resources` refers only to the NCBI resources examined by the program and should not be interpreted as proof of absence from all NCBI or INSDC resources.

### Program 4: Final cross-repository analysis

`code/04_APBT_CrossRepository_Final_Analysis.py`

Runs the final 23-check analysis using the cross-repository enrichment produced by Program 3, generates the final exposure scores, sensitivity analyses, tables, figures, and validation outputs used for the revised manuscript.

For the frozen 100,000-record manuscript cohort, the final analysis reproduces a mean exposure score of **5.22374** and **51,492 records (51.492%)** classified as High or Critical.

## Data

### BioSample retrieval dataset

The BioSample-level dataset produced by Program 1 is available in this repository:

`data/SARS_CoV2_BioSample_metadata_100000.csv`

### Cross-repository enriched dataset

The enriched cross-repository dataset used for the final manuscript analysis is archived on Zenodo because it exceeds GitHub's browser upload limit.

**Exact manuscript dataset, Version 1:**  
https://doi.org/10.5281/zenodo.22084045

**Concept DOI, representing all versions:**  
https://doi.org/10.5281/zenodo.22084044

The version-specific DOI should be used when reproducing the exact results reported in the manuscript.

## Cross-repository checks

The cross-repository enrichment is applied to:

- **P2:** BioProject association
- **P3:** Nucleotide/GenBank linkage
- **P4:** SRA linkage
- **A1:** sequencing platform or technology
- **A3:** coverage or sequencing depth
- **A4:** assembly method

Information found directly in BioSample or recovered from linked NCBI resources is treated as available. An exposure indicator is assigned only when the corresponding information remains unavailable after the queried linked resources have been examined.

## Audit terminology

The empirical analysis evaluates **metadata exposure indicators** within the APBT framework. The findings do not demonstrate deliberate manipulation, intentional poisoning, or confirmed cybersecurity compromise.

## Requirements

The programs use Python 3 and common scientific Python packages, including:

- Biopython
- pandas
- NumPy
- Matplotlib

Install required packages as needed, for example:

```bash
pip install biopython pandas numpy matplotlib
```

## NCBI Entrez usage

Programs that access NCBI Entrez require a contact email. Supply it at runtime using the relevant `--email` option or the `NCBI_EMAIL` environment variable. An NCBI API key may also be supplied where supported.

## Citation

If you use the cross-repository dataset, please cite the archived Zenodo record:

**Anjum, N. (2026). APBT Genomic Metadata Audit: Cross-Repository SARS-CoV-2 BioSample Dataset. Zenodo.**  
https://doi.org/10.5281/zenodo.22084045

The associated article should also be cited once its final bibliographic details are available.

## License

The Zenodo dataset is released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

Repository source-code licensing should be interpreted according to the license file included with this repository, if present.
