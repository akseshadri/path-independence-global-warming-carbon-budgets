# Path Independence Between Cumulative Emissions and Global Warming

Code accompanying the papers:

> 1. Seshadri, A. K. (2017). Origin of path independence between cumulative CO₂ emissions and global warming. *Climate Dynamics*, 49, 3383–3401. [doi:10.1007/s00382-016-3519-3](https://doi.org/10.1007/s00382-016-3519-3)

> 2. Seshadri, A. K. (2021). Cumulative emissions accounting of greenhouse gases due to path independence for a sufficiently rapid emissions cycle. *Climate Dynamics*. [doi:10.1007/s00382-021-05739-3](https://doi.org/10.1007/s00382-021-05739-3)

## Overview

Observations and GCMs show that global warming is approximately proportional to cumulative CO₂ emissions, regardless of the emissions pathway — a property known as "path independence". This is the foundation of cumulative emissions accounting for climate mitigation policy. These papers use a two-box energy balance model (EBM) coupled with models of atmospheric CO₂ to identify the mathematical conditions under which path independence arises.

Paper 1 identifies sufficient conditions as a pair of timescale inequalities (ζ\_r ≫ ζ\_{m₁} ≪ τ\_D), showing that slow carbon cycle processes are essential while large deep-ocean heat capacity is not. Paper 2 extends the framework to general greenhouse gases, showing that path independence holds when the emissions cycle period is comparable to or shorter than the atmospheric lifetime.

## Repository Structure

```
├── README.md
├── setup.m                                    # Run first: adds all paths to MATLAB
├── .gitignore
├── shared/                                    # Shared model functions (both papers)
│   ├── twoboxmodelco2.m                       # Two-box EBM with CO2-only forcing
│   ├── co2conc_co2only.m                      # CO2 concentration via Joos et al. (2013)
│   ├── forcing.m                              # Total radiative forcing calculator
│   └── forcingco2.m                           # CO2-only radiative forcing
├── paper1_origin_of_path_independence/        # Seshadri (2017)
│   ├── figures/
│   │   ├── PaperFig1.m                        # Fig. 1: EBM approximation verification
│   │   ├── PaperFig1_secondhalf.m             # Fig. 1 (extended scenarios)
│   │   ├── PaperFig2_SIFig1.m                 # Fig. 2 & SI Fig. 1: Airborne fraction
│   │   ├── PaperFig3.m                        # Fig. 3: θ evolution
│   │   ├── PaperFig4.m                        # Fig. 4: Timescale comparison
│   │   ├── PaperFig5_6.m                      # Figs. 5–6: Sensitivity & path independence
│   │   ├── PaperFig5_6_secondhalf.m           # Figs. 5–6 (extended scenarios)
│   │   ├── PaperFig7.m                        # Fig. 7: Effect of damping timescale
│   │   ├── PaperFig8.m                        # Fig. 8: Effect of f_long
│   │   
│   ├── si_figures/
│   │   ├── SIFig2.m                           # SI Fig. 2: τ_D spread across GCMs
│   │   ├── SIFig3.m                           # SI Fig. 3: Additional scenarios
│   │   └── SIFig4.m                           # SI Fig. 4: Airborne fraction differences
│   └── helpers/
│       ├── co2conc_co2onlyv2.m                # CO2 variant: piecewise-linear emissions
│       └── co2conc_co2onlywp2.m               # CO2 variant: multiplicative perturbation
├── paper2_cumulative_emissions_accounting/    # Seshadri (2021)
│   ├── figures/
│   │   ├── PaperFigs_1_2.m                    # Figs. 1–2: Path independence for CO2
│   │   ├── PFig3.m                            # Fig. 3: Cubic polynomial condition
│   │   ├── PFig4.m                            # Fig. 4: Tolerance contour plot
│   │   ├── PFig5.m                            # Fig. 5: |τ_M/τ_r| vs α
│   │   └── PFig6.m                            # Fig. 6: HFC-143a example
│   └── helpers/
│       ├── getconc.m                          # Concentration for single-lifetime species
│       └── twoboxmodelsingleforcer.m          # EBM for generic forcing time series
└── data/                                      # Input data files
    ├── co2.txt                                # Historical CO2 emissions
    ├── betaandnu.txt                          # EBM parameters from CMIP5
    ├── heatcapacity.txt                       # Heat capacity estimates
    ├── histF.txt                              # Historical radiative forcing
    ├── hfchist.xlsx                           # Historical HFC emissions
    └── rcp*proj.xlsx                          # RCP scenario projections
```

## Requirements

- **MATLAB** (tested with R2016b and later)
- No additional toolboxes required

## Usage

1. Open MATLAB and run `setup` from the repository root to add all function paths:
   ```matlab
   cd('/path/to/path-independence-ebm')
   setup
   ```

2. Navigate to a figure directory and run any script:
   ```matlab
   cd paper1_origin_of_path_independence/figures
   PaperFig1
   ```

   Data files are placed in each figure directory so that `textread`/`xlsread` calls resolve in the current working directory. Functions are resolved via the MATLAB path configured by `setup.m`.

## Citation

If you use this code, please cite the relevant paper(s):

```bibtex
@article{Seshadri2017,
  title   = {Origin of path independence between cumulative {CO}$_2$ emissions and global warming},
  author  = {Seshadri, Ashwin K.},
  journal = {Climate Dynamics},
  volume  = {49},
  pages   = {3383--3401},
  year    = {2017},
  doi     = {10.1007/s00382-016-3519-3}
}

@article{Seshadri2021,
  title   = {Cumulative emissions accounting of greenhouse gases due to path independence for a sufficiently rapid emissions cycle},
  author  = {Seshadri, Ashwin K.},
  journal = {Climate Dynamics},
  year    = {2021},
  doi     = {10.1007/s00382-021-05739-3}
}
```

## Contact

Ashwin K. Seshadri — [ashwins@iisc.ac.in](mailto:ashwins@iisc.ac.in)
Indian Institute of Science, Bangalore 560012, India
