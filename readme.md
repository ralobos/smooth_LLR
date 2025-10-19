# Smooth Local Low-Rank (LLR) Reconstruction Software v1.0

## Overview

This MATLAB software reproduces the reconstruction experiments presented in [1]. It performs reconstruction of retrospectively undersampled dynamic MRI k-space data using a **smooth Huber-based local low-rank regularizer**. The key innovation is the use of smooth regularizers that enable standard optimization algorithms (such as nonlinear conjugate gradient) to solve the inverse problem efficiently.

## Contents

### Main Scripts
- **`example_smooth_LLR.m`** - Main reconstruction script for dynamic MRI undersampled data
- **`example_smooth_LLR_fast_step_size.m`** - Reconstruction using a heuristic fast step-size selection strategy


## References

**[1]** R. A. Lobos, J. Salazar Cavazos, R. R. Nadakuditi, J. A. Fessler.  
*Smooth optimization algorithms for global and locally low-rank regularizers*  
arXiv:2505.06073, May 2025.


