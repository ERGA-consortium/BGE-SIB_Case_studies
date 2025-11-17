# A new reference genome for the conservation of coastal species in the face of anthropogenic change. _Dysdera curvisetae_ Genome Assembly and Annotation Pipeline

This repository contains a complete pipeline for genome assembly, scaffolding, quality assessment, repeat masking, gene prediction, and functional annotation for _Dysdera curvisetae_ or related organisms.  
The pipeline is designed to be adapted to be run on HPC environments and uses PacBio HiFi reads, Hi-C data, RNA-seq, and protein evidence.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Pipeline Overview](#pipeline-overview)
---

## Prerequisites

- **Tools**:
  - `hifiasm 0.24.0-r702`
  - `bwa 0.7.17-r1188`, `samtools 1.19.2`, `HapHiC 1.0.7`
  - `BUSCO 5.8.0 arachnida_odb10`
  - `QUAST v5.0.2`
  - `meryl 1.3` and `Merqury v1.3`
  - `RepeatModeler 2.0.3` & `RepeatMasker 4.1.7-p1`
  - `BRAKER3 3.0.8`, `TSEBRA v1.1.2.5`
  - `BLAST 2.13.0+`
  - `emapper-2.1.3` (eggNOG)
  - `InterProScan 5.72-103.0`


- **Data**:
  - PacBio HiFi reads
  - Hi-C reads
  - RNA-seq datasets
  - Curated protain database for annotation
  - Protein reference databases (Arthropoda DB, SwissProt)

---

## Pipeline Overview

The pipeline follows these steps, with scripts provided for each stage:

1. **HiFi + Hi-C Assembly**  
   Assemble PacBio HiFi reads and scaffold with Hi-C data.  
   Script: `1-hifiasm_with_hic.sh`

2. **Hi-C Read Mapping & Filtering**  
   Align Hi-C reads to the assembly, filter BAM files and scaffolding.
   Scripts:  
   - `2.1-alignment.sh` (map reads)  
   - `2.2-filtering.sh` (filter alignments)  
   - `2.3-HapHiC_running.sh` (run HapHiC scaffolding)

3. **Assembly Quality Assessment – QUAST**  
   Evaluate basic assembly statistics using QUAST.  
   Script: `4-quast.sh`

4. **Assembly Quality Assessment – BUSCO**  
   Evaluate assembly completeness using BUSCO.  
   Script: `3-BUSCO.sh`  
   > Note: BUSCO can also be run after gene prediction to assess annotation quality.

5. **K-mer Quality Assessment**  
   Assess assembly completeness and error profiles using Merqury.  
   Script: `5-merqury.sh`

6. **Repeat Identification and Masking**  
   Detect and mask repetitive elements using RepeatModeler and RepeatMasker.  
   Script: `6-RepeatModeler_Masker.sh`

7. **Gene Prediction**  
   Predict genes using BRAKER3, refine models with TSEBRA, and extract longest isoforms.  
   Scripts:  
   - `7.1-BRAKER3.sh` (initial gene prediction)  
   - `7.2-BRAKER3_extra_step.sh` (TSEBRA refinement)

8. **Functional Annotation**  
   Annotate predicted proteins using BLASTP, eggNOG, and InterProScan.  
   Scripts:  
   - `8.1-Functional_BLASTP_ArthropodaDB.sh`  
   - `8.1-Functional_BLASTP_SwissProt.sh`  
   - `8.1-Functional_eggnog.sh`  
   - `8.1-Functional_InterProScan.sh`

---

## Notes

- Each script is self-contained and can be run independently once the required inputs from previous steps are ready.  
- Adjust CPU/thread numbers and file paths according to your HPC environment.  
- The scripts are numbered to reflect the recommended order of execution, but you can run some steps in parallel if resources allow.  

---

## Author

Sara Guirao-Rico, University of Barcelona, sguirao@ub.edu  
Nuria Macías-Hernández, University of La Laguna, nemacias@ull.edu.es  
Julio Rozas, University of Barcelona, jrozas@ub.edu  
Silvia Garcia-Juan, University of Barcelona, silvia.garciaj@ub.edu  

# Acknowledgements
This work was supported by the SIB Swiss Institute of Bioinformatics under the Biodiversity Genomics Europe project.


