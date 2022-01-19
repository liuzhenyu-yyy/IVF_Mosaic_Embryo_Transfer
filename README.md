# IVF Mosaic Embryo Transfer

 Custom code for IVF project.

We applied single-cell multiomics sequencing for seven infants with blastula chromosomal mosaicism detected by TE biopsy.

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

