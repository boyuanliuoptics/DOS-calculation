# Density-of-states calculation in photonic crystals

This is a Python implementation by Bradley Dice (@bdice, University of Michigan) of Scheme and MATLAB code originally written by Boyuan Liu.
Please see the article ['Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals'](https://doi.org/10.1088/2040-8986/aaae52) for more information.

## Requirements

This code requires the following libraries.

1. ase for generating k-points.
2. MPB (with Python bindings) for computing photonic bands.
3. Matplotlib for plotting data.
4. NumPy for numerical computations.
5. tqdm for showing progress bars.

I recommend to install all of the libraries from conda-forge:
```bash
conda install -c conda-forge ase mpb numpy tqdm
```
