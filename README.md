# MEDIOLA-WOO

Analysis code and processed data supporting the MEDIOLA window-of-opportunity (WOO)
trial of olaparib and durvalumab in triple-negative or estrogen receptor (ER)-low
breast cancer (ClinicalTrials.gov [NCT03594396](https://clinicaltrials.gov/study/NCT03594396)).

## Study overview

MEDIOLA-WOO was a single-arm window-of-opportunity trial in 54 patients
with triple-negative or ER-low breast cancer. Patients received olaparib
continuously for 4 weeks, with durvalumab added at week 2, before standard
neoadjuvant chemotherapy. Paired tumor biopsies were collected at baseline,
on-treatment week 2 (OT1), and week 4 (OT2), together with FDG-PET/MRI and CT
imaging. Using whole-exome sequencing, bulk RNA-seq, AI-based digital pathology,
and multiplex immunohistochemistry, we characterized the early immunomodulatory
effects of PARP inhibition and their relationship to treatment response, to
inform the rationale and sequencing of PARP-inhibitor plus immune-checkpoint
blockade combinations.

## Repository structure

```
MEDIOLA-WOO/
├── scripts/     R scripts that generate each main and supplementary figure and table
├── data/        processed, study-generated input files (see data/README.md)
│   └── README.md   describes every file, its source, and which script consumes it
├── README.md
└── LICENSE
```

## Reproducing the figures and tables

1. Download the sequencing-derived inputs listed under **Data availability** into `data/`.
2. Open a script in `scripts/`, set `data_dir` to your local `data/` folder, and run it in
   R (≥ 4.5). Required packages are listed at the top of each script.
3. Outputs (PDFs and result tables) are written to the working directory.

`data/README.md` lists, for each input file, whether it is downloaded from a public
archive, rebuilt from a public reference, or provided here, and which script uses it.

## Data availability

- **Processed RNA-seq** — Gene Expression Omnibus, accession
  [`GSE341672`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341672).
- **Raw RNA-seq and whole-exome sequencing** — European Genome-phenome Archive,
  study accession `EGAS50000002019`, comprising datasets `EGAD50000002906`
  (RNA sequencing) and `EGAD50000002907` (whole-exome sequencing). These data
  contain potentially identifying germline genetic information and are therefore
  under controlled access; requests are administered by the MEDIOLA WOO Data
  Access Committee (`EGAC50000001066`). Requests must be submitted through the
  EGA website.
- **Source Data** for all figures are provided with the published article.

## Code availability

All custom analysis and figure-generation code is contained in this repository.
Upstream read alignment, variant calling, and expression quantification were performed
with standard tools at the versions and parameters given in the Methods.

## Citation

If you use this code or data, please cite the archived release:
https://doi.org/10.5281/zenodo.21802544

## Contact

Questions about the code can be raised as a GitHub issue or sent to Jiwon Koh
(jwskoh@snu.ac.kr), Department of Pathology, Seoul National University Hospital.
Correspondence regarding the study should be addressed to Seock-Ah Im.

## License

Code in this repository is released under the MIT License (see `LICENSE`).
