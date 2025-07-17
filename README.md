# Unequal Randomization in Clinical Trials

This repository contains R scripts for simulating and analyzing different randomization schemes in two-arm clinical trials. The simulations focus on evaluating **key performance metrics** such as expected mean response, Type-I error rate, statistical power, and patient allocation across various trial designs.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Project Context](#project-context)
- [Randomization Schemes Implemented](#randomization-schemes-implemented)
- [Key Metrics Analyzed](#key-metrics-analyzed)
- [Files in this Repository](#files-in-this-repository)
- [Reproducibility](#reproducibility)
- [Dependencies](#dependencies)

---

## Project Overview

This project provides a simulation framework to compare the performance of various **randomization strategies** for two-arm clinical trials. It lets researchers investigate how different allocation rules impact trial outcomes under diverse data distributions and parameter settings.

---

## Project Context

The simulations in this repository are designed to generate data for specific tables and figures in a research study:

* **Table 1**: Results from simulations using **binary outcomes** (i.e., `dist="bern"` in `sim_ER`, `sim_FR`, `sim_ERADE`) are intended for Table 1.
* **Table 2**: Results from simulations using **normal outcomes** (i.e., `dist="norm"` in `sim_ER`, `sim_FR`, `sim_ERADE`), as shown in the examples below, are intended for Table 2.
* **Figure 1**: The **allocation ratios** (specifically, the mean proportion of patients allocated to the experimental arm, corresponding to column 6 of the simulation output) from all simulation types are relevant for illustrating patient distribution and are intended for Figure 1.

---

## Randomization Schemes Implemented

The simulations cover the following randomization methods:

1.  **Equal Randomization (ER)**: Patients are assigned to treatment arms with equal probability (e.g., 50/50).
2.  **Fixed Randomization (FR)**: Patients are assigned to treatment arms with a pre-specified fixed probability (e.g., 70% to arm 1, 30% to arm 2).
3.  **ERADE (Efficient Response Adaptive Design)**: An adaptive randomization method that adjusts allocation probabilities based on observed patient responses during the trial. This includes specific implementations like "RSHIR_Z0", "RSHIR_Z1", "Neyman_Z0", and "Neyman_Z1," which adapt based on different variance and mean estimators.

---

## Key Metrics Analyzed

For each simulation, we compute and evaluate the following metrics:

* **Expected Mean Response (EMR)**: The overall average response observed across all patients in a trial, averaged over many simulations.
* **Type-I Error Rate**: The probability of incorrectly rejecting the null hypothesis (i.e., concluding there's a difference between arms when there is none). We assess this when there's no true difference between the experimental and control arms.
* **Power**: The probability of correctly rejecting the null hypothesis (i.e., detecting a true difference between arms when one exists). We assess this when a true difference is present.
* **Allocation towards Experimental Arm**: The average proportion of patients assigned to the experimental treatment arm.
* **Proportion in Theoretically Superior Arm**: For adaptive designs, this indicates how often patients were allocated to the arm that's truly better.

---

## Files in this Repository

* `Data_Generator.R`: Contains functions to simulate patient data based on specified distributions (e.g., Bernoulli, Normal) and parameters. It also likely includes the `superior()` function to determine which arm is theoretically superior based on the true parameters.
* `WaldTest.R`: Provides functions (`wald.test`, `wald.test.binary`) for performing Wald tests to compare outcomes between two arms.
* `sim_ER.R`: The main simulation function for Equal Randomization. It calls `two_arm_ER`.
* `two_arm_ER.R`: Implements a single trial simulation for Equal Randomization.
* `sim_FR.R`: The main simulation function for Fixed Randomization. It calls `two_arm_FR`.
* `two_arm_FR.R`: Implements a single trial simulation for Fixed Randomization.
* `sim_ERADE.R`: The main simulation function for ERADE (Adaptive Randomization). It calls `two_arm_ERADE`.
* `two_arm_ERADE.R`: Implements a single trial simulation for ERADE, including the logic for dynamic patient allocation.
* **Helper Functions**: The `two_arm_ERADE` script relies on helper functions like `ERADE()` (to calculate ERADE probabilities) and `find_root()` (for RSHIR_Z0 calculations), which are included in the same file or implicitly sourced.

---

## Reproducibility

All results presented in the associated paper can be generated using the respective simulation files and R code provided within this repository. This ensures full **transparency and reproducibility** of the research findings.

---

## Dependencies

This project requires the following R packages:

* `BSDA`
* `rankFD`
* `lawstat`

You can install them using:

```R
install.packages(c("BSDA", "rankFD", "lawstat"))
