# Origin of Path Independence Between Cumulative CO₂ Emissions and Global Warming

Code accompanying the paper:

> Seshadri, A. K. (2017). Origin of path independence between cumulative CO₂ emissions and global warming. *Climate Dynamics*, 49, 3383–3401. [doi:10.1007/s00382-016-3519-3](https://doi.org/10.1007/s00382-016-3519-3)

## Overview

This paper identifies sufficient conditions for path independence using a closed-form expression for global warming in a two-box EBM. Path independence arises from weak inequality constraints: ζ\_r(t) ≫ ζ\_{m₁}(t) ≪ τ\_D, where ζ\_r is the timescale for airborne fraction changes, ζ\_{m₁} is the cumulative emissions timescale, and τ\_D is the damping timescale. Slow carbon cycle time-constants are essential; large deep-ocean heat capacity is not.

## Figure Descriptions

| Script | Paper Figure | Description |
|--------|-------------|-------------|
| `PaperFig1.m` | Fig. 1 | Verification of EBM approximation (Eq. 1) versus numerical integration |
| `PaperFig1_secondhalf.m` | Fig. 1 (cont.) | Extended emissions scenarios for EBM verification |
| `PaperFig2_SIFig1.m` | Fig. 2 & SI Fig. 1 | Airborne fraction dynamics: contributions h₁, h₂, h₃ from Eq. (3)–(4) |
| `PaperFig3.m` | Fig. 3 | Evolution of θ (ratio of warming to atmospheric CO₂ increase) |
| `PaperFig4.m` | Fig. 4 | Comparison of timescales ζ\_r, ζ\_{m₁}, τ\_D for path independence |
| `PaperFig5_6.m` | Figs. 5–6 | Sensitivity dT\_s/dm₁ and global warming vs cumulative CO₂ |
| `PaperFig5_6_secondhalf.m` | Figs. 5–6 (cont.) | Extended emissions scenarios |
| `PaperFig7.m` | Fig. 7 | Effect of damping timescale τ\_D — path independence is robust |
| `PaperFig8.m` | Fig. 8 | Effect of f\_long — slow carbon cycle is essential |
| `SIFig2.m` | SI Fig. 2 | Spread of τ\_D across 16 CMIP5 GCMs |
| `SIFig3.m` | SI Fig. 3 | Additional scenarios (late-peaking, second phase) |
| `SIFig4.m` | SI Fig. 4 | Airborne fraction differences across scenarios |

## Usage

```matlab
cd('/path/to/path-independence-ebm')
setup
cd paper1_origin_of_path_independence/figures
PaperFig1
```
