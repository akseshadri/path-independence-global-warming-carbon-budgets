# Cumulative Emissions Accounting of Greenhouse Gases for Sufficiently Rapid Emissions Cycle

Code accompanying the paper:

> Seshadri, A. K. (2021). Cumulative emissions accounting of greenhouse gases due to path independence for a sufficiently rapid emissions cycle. *Climate Dynamics*. [doi:10.1007/s00382-021-05739-3](https://doi.org/10.1007/s00382-021-05739-3)

## Overview

This paper extends the path independence framework beyond CO₂ to greenhouse gases with a single atmospheric lifetime τ. Path independence depends on the ratio α = T/τ between the emissions cycle period and the atmospheric lifetime, being valid when T ≲ τ. This makes cumulative emissions accounting potentially relevant to HFCs and other GHGs with lifetimes of several decades whose emissions have recently begun.

## Figure Descriptions

| Script | Paper Figure | Description |
|--------|-------------|-------------|
| `PaperFigs_1_2.m` | Figs. 1–2 | Path independence conditions for CO₂: directional derivatives and timescale bounds |
| `PFig3.m` | Fig. 3 | Cubic polynomial g(y) determining path independence condition (Eq. 32) |
| `PFig4.m` | Fig. 4 | Contour plot of tolerance θ versus α = T/τ and rescaled time x = t/T |
| `PFig5.m` | Fig. 5 | Ratio \|τ\_M/τ\_r\| at emissions peak versus α, for different shape functions |
| `PFig6.m` | Fig. 6 | HFC-143a example: slow vs rapid mitigation scenarios |

## Usage

```matlab
cd('/path/to/path-independence-ebm')
setup
cd paper2_cumulative_emissions_accounting/figures
PFig6
```
