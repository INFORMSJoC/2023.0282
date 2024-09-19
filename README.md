![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)

# Reproducible feature selection for high-dimensional measurement error models

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](https://github.com/INFORMSJoC/2023.0282/blob/master/LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper
[Reproducible feature selection for high-dimensional measurement error models](https://doi.org/10.1287/ijoc.2023.0282)
by X. Zhou, Y. Li, Z. Zheng, J. Wu and J. Zhang.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using the following DOIs.

[https://doi.org/10.1287/ijoc.2023.0282](https://doi.org/10.1287/ijoc.2023.0282)

[https://doi.org/10.1287/ijoc.2023.0282.cd](https://doi.org/10.1287/ijoc.2023.0282.cd)

Below is the BibTex for citing this version of the code.
```latex
@misc{zhou2024reproducible,
  author =        {X. Zhou, Y. Li, Z. Zheng, J. Wu and J. Zhang},
  publisher =     {INFORMS Journal on Computing},
  title =         {Reproducible feature selection for high-dimensional measurement error models},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0282.cd},
  url =           {https://github.com/INFORMSJoC/2023.0282},
}  
```

## Description

This directory contains the code for the *double projection knockoff filter (DP-knockoff)* algorithm.

This project contains four folders: `data`, `results`, `src`, `scripts`.
- `data`: include datasets used in the paper.
- `results`: include the experimental results.
- `src`: include the source codes.
- `scripts`: include codes to replicate the experiments in the paper.

## Replicating

To get the Figures and Tables in the paper, please run the codes in the `scripts` folder. 
- `scripts/Tables 1-2`: to get Tables 1-2 of the paper
- `scripts/Tables 3-4`: to get Tables 3-4 of the paper
- `scripts/Figures 1-3`: to get Figures 1-3 of the paper
- `scripts/Tables 5-6`: to get Tables 5-6 of the paper
- `scripts/Tables 7-8`: to get Tables 7-8 of the paper

## Data

In Section 6.1, the complete dataset is available at: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MEXP-1618?query=E-MEXP-1618. The `data` only contains one file 'Pasient2.CEL' due to the large file size.

## Ongoing Development

This code is being developed on an on-going basis at the author's [GitHub site](https://github.com/xinaut/2023.0282).
