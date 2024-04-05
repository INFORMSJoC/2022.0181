# 2022.0181

### ![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)

# Robust Parallel Pursuit for Large-Scale Association Network Learning
This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](https://github.com/INFORMSJoC/2019.0000/blob/master/LICENSE).
## Cite
To cite the contents of this repository, please cite both the paper and this repo, using the following DOIs.

[https://doi.org/10.1287/ijoc.2022.0181](https://doi.org/10.1287/ijoc.2022.0181)

[https://doi.org/10.1287/ijoc.2022.0181.cd](https://doi.org/10.1287/ijoc.2022.0181.cd)

Below is the BibTex for citing this version of the code.
```latex
@article{Li2024IJOC,
  author =        {W. Li, X. Zhou, R. Dong, and Z. Zheng},
  publisher =     {INFORMS Journal on Computing},
  title =         {Robust Parallel Pursuit for Large-Scale Association Network Learning, v2022.0181},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0181.cd},
  url =           {https://github.com/INFORMSJoC/2022.0181},
}  
```
## Description
Sparse reduced-rank regression is an important tool to uncover the large-scale response-predictor association network,  as exemplified by modern applications such as the diffusion networks, and recommendation systems. However, the association networks recovered by existing methods are either sensitive to outliers or not scalable under the big data setup. In this paper, we propose a new statistical learning method called robust parallel pursuit (ROP) for joint estimation and outlier detection in large-scale response-predictor association network analysis.  The proposed method is scalable in that it transforms the original large-scale network learning problem into a set of sparse unit-rank estimations via factor analysis,  thus facilitating an effective parallel pursuit algorithm.  

This project contains four folders: `data`, `results`, `src`, `scripts`.
- `data`: include two datasets used in the paper.
- `results`: include the experimental results.
- `src`: include the source codes.
- `scripts`: include codes to directly replicate the experiments in the paper.

## Replicating
To get the Figures and Tables in the paper, please run the codes in the `scripts` folder. 
- `scripts/Table1 Figures1-5`: to get Table 1 and Figures 1-5 of the paper
- `scripts/Table2 Figure6`: to get Table 2 and Figure 6 of the paper
- `scripts/Table3`: to get Table 3 of the paper
- `scripts/Table4`: to get Table 4 of the paper

Specifically, for example, in `scripts/Table1 Figures1-5`,
- Figure1.R: run Figure1.R in this file to get Figure 1 of the paper
- Figure2.R: run Figure2.R in this file to get Figure 2 of the paper
- Figure3.R: run Figure3.R in this file to get Figure 3 of the paper
- Figure4.R: run Figure4.R in this file to get Figure 4 of the paper
- Figure5.R: run Figure5.R in this file to get Figure 5 of the paper
- Table1.R: run Table1.R in this file to get Table 1 of the paper
  

See the README.md file in each folder for a detailed description.

## Ongoing Development
This code is being developed on an on-going basis at the author's [Github site](https://github.com/zeminzhe/2022.0181).