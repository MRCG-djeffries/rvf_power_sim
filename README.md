# RVF Vaccine Power Simulations

**Exact Binomial Case-Split · Event-Driven · Group-Sequential Design**

This repository contains the R code used to generate the power and sample size simulations described in **statsbits3.docx**, which supports the statistical methods and sample size justification for the RVF vaccine trial protocol.

The simulations are event-driven, use an exact discrete binomial case-split test, and incorporate group-sequential monitoring, latent clustering of infection risk, and baseline seropositivity. The code was executed on a high-performance computing (HPC) system using SLURM, but can also be run locally for debugging and testing.

---

## Overview

The simulations are designed to determine:

- the number of seronegative endpoint events required to achieve target power across a range of vaccine efficacies (VE), and  
- the corresponding seronegative and total enrollment (expressed in person-years) required under different assumptions about baseline seroprevalence and transmission heterogeneity.

---

## Repository structure

- **run_test_power9HPC.R**  
  Main simulation driver used for all production runs.

- **slurm_run_test_power9HPC_binsplit.sh**  
  SLURM submission script used to run the simulations on the HPC.

- **plots_out1.R**  
  Code used to generate the figures included in *statsbits3.docx*.

- **out5/**  
  Poisson-based event-count estimates used to initialise and accelerate the event-number search for the exact binomial case-split design.

- **out11/**  
  Final simulation outputs from the exact binomial case-split analysis.

The directory **out11** is included so that all figures and results in the document can be reproduced without rerunning the full HPC simulations.

---

## Parameter grid

The simulation parameter grid contains **242 rows**:

  ## -------- define parameter grid (2×11×11 = 242) --------
  wave_pct <- c(0.1,1)           
  
  hr_true  <- rev(seq(0.05,0.25,length.out=11))#c(0.25, 0.20, 0.15, 0.10, 0.05)
  
  sero_pct <- seq(10,50,length.out=11)#c(10, 20, 30, 40, 50)          
  
  grid <- expand.grid(
    wave = wave_pct,
    hr   = hr_true,
    sero = sero_pct,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

Each SLURM array task evaluates one row of this grid.
where each row corresponds to a unique combination of:

Attack rate

Vaccine efficacy

Baseline seroprevalence

Transmission heterogeneity assumptions

---

## Running the simulations on the HPC

Simulations are run using SLURM via:

```bash
slurm_run_test_power9HPC_binsplit.sh
```
## Running the Code Locally (Debugging)

The same code can be run directly from an R session on any machine for debugging or testing.

Example:

run_test_power5(242, 100, "out11", "binsplit", "out5")

Arguments

242
Row of the parameter grid (1–242)

100
Number of simulation replicates used to estimate power

"out11"
Output directory for binomial case-split results

"binsplit"
Analysis type (alternative is "poisson")

"out5"
Directory containing Poisson-based event-number estimates

The Poisson-based estimates are used to seed the search for required event numbers in the exact binomial case-split simulations. A simple binary search is not possible because power as a function of events exhibits a saw-tooth pattern due to discreteness of the binomial test.

## Serostatus–Risk Correlation

The correlation between baseline serostatus and infection risk is controlled by the parameter, b_sero

This parameter governs the strength of enrichment of seropositivity in high-risk clusters.

Note:
b_sero was introduced late in development and is currently hardcoded in
sim_mult_wave8_plus() at approximately line 539 of run_test_power9HPC.R.

## Software Requirements

R version:
Flexible; simulations were run on the HPC using R 4.3.2 and desktop testing used R 4.3.3

Required R packages:

gsDesign (version 3.8.0)

No other non-base packages are required.

## Reproducibility

All figures and tables in statsbits3.docx can be reproduced using:

The included out5 and out11 result directories

The plotting script plots_out1.R

Full reruns of the simulations are computationally intensive and were performed on an HPC environment.
