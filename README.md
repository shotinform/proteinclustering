# Protein Clustering: PTA–CVD–T2D Network-Medicine Replication

This repository contains code, data, and analysis notebooks for replicating the study:

> *Network-medicine approach for the identification of genetic association of parathyroid adenoma with cardiovascular disease and type-2 diabetes*  
> (Nikhat Imam et al., Briefings in Functional Genomics, 2023, 22, pp. 250–262)

The project was completed as part of a seminar at the Faculty of Mathematics, University of Belgrade (2023/24)

**Authors:** Marko Nikitović, Boško Andrić  
**Mentor:** Prof. Dr. Nenad Mitić  

---

## Objectives

- Collect proteins associated with **parathyroid adenoma (PTA)**, **cardiovascular disease (CVD)**, and **type-2 diabetes (T2D)**.  
- Build a **protein–protein interaction (PPI) network** using BioGRID and analyze its topological properties.  
- Detect **functional modules** using the MCODE algorithm.  
- Perform **gene set enrichment** (GO terms and Reactome pathways).  
- Explore **drug–target and disease–drug interactions** using COREMINE.

## Results (Highlights)

- **Common genes:** 33 overlapping proteins (per original study; our DisGeNET run gave 29).  
- **PPI network:** ~7,075 nodes, ~192,650 edges after filtering.  
- **MCODE modules:** 11 submodules with hub proteins (e.g., TP53, ESR1, CTNNB1, EGFR).  
- **GO/Pathway enrichment:** Significant BP/MF terms and Reactome pathways; consistent with the reference study.  
- **Drug mappings:** Identified Target-Associated Drugs (TAD) and Disease-Associated Drugs (DAD); overlap visualized with Venns and Cytoscape networks.  

---

## References

- [DisGeNET](https://www.disgenet.org/) – disease gene sets  
- [BioGRID](https://thebiogrid.org/) – protein–protein interactions  
- [IntAct](https://www.ebi.ac.uk/intact/home) – molecular interaction database  
- [Cytoscape](https://cytoscape.org/) – network visualization  
- [ClusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) – GO enrichment  
- [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) – pathway enrichment  
- [COREMINE](https://www.coremine.com/) – drug–gene/disease data  
