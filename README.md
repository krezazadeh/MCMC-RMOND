# MCMC-RMOND
**Relativistic MOND Parameter Estimation via Parallel MCMC**

---

## Overview
**MCMC-RMOND** is a Fortran 90 program that performs **Markov Chain Monte Carlo (MCMC)** sampling using the **Metropolis–Hastings** algorithm to estimate the parameters of the **Relativistic MOND (Modified Newtonian Dynamics)** model from galactic rotation curve data.

The program compares the fits of:
- **Relativistic MOND**,
- **Newtonian dynamics**, and
- **Dark Matter** models

to the observed **outer stellar rotation curves** of galaxies, such as NGC 3198.

Parallelization is implemented using **MPI**, allowing multiple chains to run simultaneously.
Convergence is evaluated using the **Gelman–Rubin** \( \hat{R} \) diagnostic.

---

## Features
- Parallel **MCMC** sampling using **MPI**
- **Metropolis–Hastings** acceptance criterion
- Automatic convergence monitoring with **Gelman–Rubin**
- Flexible **priors**, **jump sizes**, and **initial parameters**
- Comparison of **Relativistic MOND**, **Newtonian**, and **Dark Matter** fits
- Modular and well-documented **Fortran 90** structure

---
