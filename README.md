# IVF Mosaic Embryo Transfer

 Custom code for IVF project.

We applied single-cell multiomics sequencing for seven infants with blastula chromosomal mosaicism detected by TE biopsy.

For more information, please check our [publication](https://www.sciencedirect.com/science/article/pii/S1672022922000882) at Genomics, Proteomics & Bioinformatics:
> Gao Y, Zhang J, Liu Z, et al. Single-cell Sequencing Reveals Clearance of Blastula Chromosomal Mosaicism in In Vitro Fertilization Babies [published online ahead of print, 2022 Aug 6]. *Genomics Proteomics Bioinformatics*. 2022;S1672-0229(22)00088-2. doi:10.1016/j.gpb.2022.07.004

Structure of the  Repositories:

```
IVF_Mosac_Embryo_Transfer/
├── README.md
├── pipeline_TE_PGT_A (shell script to process CNV-Seq of TE biopsy)
│   └── CNVSeq.one.sh
├── pipeline_scMALBAC (shell script to process single-cell WGS of PBMC)
│   └── RUN.multiplexed_MALBAC.sh
├── pipeline_scRNA (shell script to process single-cell RNA-Seq of PBMC)
│   ├── S00.Get_Params.sh
│   ├── S01.scRNA_UMI_QC.sh
│   ├── S02.scRNA_UMI_STAR.sh
│   ├── S03.scRNA_UMI_UMI.sh
│   ├── S04.Check_Test.R
│   ├── scRNA_Normalization.R
│   └── scRNA_UMI.sh
├── scDNA.R (Downstream analysis of single-cell CNV)
└── scRNA.R (Downstream analysis of single-cell RNA-Seq)

3 directories, 12 files
```

