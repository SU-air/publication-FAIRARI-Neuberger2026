# Importance of hydrated aerosol particles for aerosol–fog relationships in the Italian Po Valley

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5194%2Fegusphere--2025--5419-green)](https://doi.org/10.5194/egusphere-2025-5419)
[![Status: Preprint](https://img.shields.io/badge/Status-Preprint-orange)](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-5419/)

> **Neuberger et al. (2025)** — *Preprint under review at Atmospheric Chemistry and Physics*
> 🔗 https://doi.org/10.5194/egusphere-2025-5419

---

## Overview

This repository contains the analysis code and figure scripts supporting the paper. The study quantifies the impact of **hydrated (non-activated) aerosol particles** on fog microphysical properties and visibility in the **Po Valley**, one of Europe's most polluted regions. It draws on detailed aerosol–fog observations from the **FAIRARI campaign** (San Pietro Capofiume, Italy, November 2021 – May 2022), combined with κ-Köhler hygroscopic growth theory and Large Eddy Simulation (LES) modelling with the MIMICA model.

Key findings:
- Hydrated aerosol particles can reach several micrometres in diameter, substantially biasing fog microphysical retrievals
- Excluding hydrated aerosols leads to an 81 % increase in effective diameter and an 87 % decrease in cloud droplet number concentration
- Hydrated particles contribute on average 21 % to liquid water content and account for 36 % of sub-kilometre visibility events without droplet activation
- Fog prediction in LES is strongly sensitive to the largest dry aerosol particles present

---

## Repository Structure

```
/
├── README.md
├── LICENSE
│
├── data_processing/        # Core analysis scripts
│   ├── kappa_kohler/       # κ-Köhler hygroscopic growth calculations
│   ├── microphysics/       # Fog microphysical property processing (ED, CDNC, LWC)
│   ├── visibility/         # Visibility calculations and hydrated aerosol contributions
│   └── GFAS/              # GFAS instrument data reading and processing
│
├── LES/                    # Large Eddy Simulation (MIMICA) input/output processing
│   └── ...
│
└── figures/                # Scripts to reproduce all figures in the paper
    ├── figure_01/
    ├── figure_02/
    └── ...
```

---

## Data

The observational dataset from the FAIRARI campaign is published separately:

> Neuberger et al. — *The fog and aerosol interaction research Italy (FAIRARI) campaign, November 2021 to May 2022*
> 🔗 https://bolin.su.se/data/fairari-2021-2022-3

The LES model code (MIMICA v5) used in this study is available at:
> Savre et al. — *MIMICA LES model ver 5*
> 🔗 https://bitbucket.org/matthiasbrakebusch/mimicav5/src/master

---

## Getting Started

```bash
git clone https://github.com/SU-air/publication-FAIRARI-Neuberger2026.git
cd publication-FAIRARI-Neuberger2026
```

---

## Citation

If you use this code in your work, please cite:

> Neuberger et al. (2025): *Importance of hydrated aerosol particles for aerosol–fog relationships in the Italian Po Valley*, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-5419

---

## License

This repository is licensed under the [MIT License](LICENSE).

---

## Contact

*Paul Zieger, Stockholm University, Dept of Env Science, paul.zieger@aces.su.se*
