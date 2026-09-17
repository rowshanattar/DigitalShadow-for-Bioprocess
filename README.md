# EKF-Based Digital Twin for Bioreactor State Estimation 

A Python implementation of an Extended Kalman Filter (EKF) Digital Shadow
for baker's yeast batch cultivation, developed as part of a Master's thesis
on predictive maintenance and root cause analysis in bioprocesses.

---

## Overview

This repository implements and extends the method described in:

> Yousefi-Darani, A., Paquet-Durand, O., & Hitzmann, B. (2020).
> *The Kalman Filter for the Supervision of Cultivation Processes.*
> In: Herwig, C., Pörtner, R., Möller, J. (eds) **Digital Twins.**
> Advances in Biochemical Engineering/Biotechnology, vol 177.
> Springer, Cham.
> https://doi.org/10.1007/10_2020_145

The Digital Twin combines a nonlinear biological process model (Monod
kinetics) with a Continuous-Discrete Extended Kalman Filter to
continuously estimate unmeasured state variables — biomass, glucose,
and maximum specific growth rates — from sparse ethanol measurements
arriving every 5 minutes.

---

## Repository structure

```
ekf_bioreactor/
│
├── ekf_bioreactor.ipynb        # v1.0 — baseline implementation
│                                 Reproduction of Yousefi-Darani et al.
│
├── README.md                   # this file
└── LICENSE                     # MIT License
```

New notebooks will be added for each research hypothesis:

```
ekf_bioreactor_H4a_noise.ipynb       # H4a: sensor noise robustness
ekf_bioreactor_H4b_drift.ipynb       # H4b: sensor drift robustness
ekf_bioreactor_H2_sensors.ipynb      # H2:  additional measurements
ekf_bioreactor_H1_observability.ipynb # H1: observability limits
ekf_bioreactor_H3_fault_id.ipynb     # H3: fault identification
```

---

## The Digital Twin

The EKF estimates 5 state variables simultaneously:

| Variable | Description | Unit | Measured? |
|---|---|---|---|
| X | Biomass concentration | g/L | No — estimated |
| G | Glucose concentration | g/L | No — estimated |
| E | Ethanol concentration | g/L | Yes — every 5 min |
| μmax,G | Max. growth rate on glucose | h⁻¹ | No — estimated |
| μmax,E | Max. growth rate on ethanol | h⁻¹ | No — estimated |

The process model is based on Monod kinetics (Equation 19,
Yousefi-Darani et al., 2020). The Jacobian matrix required for
covariance propagation is computed analytically using symbolic
differentiation via `sympy`.

---

## Research hypotheses

This project investigates four hypotheses beyond the baseline:

**H1 — Observability limit:**
How many unmeasured state variables can be added before the EKF
loses the ability to reliably estimate them from a single measurement?

**H2 — Additional measurements:**
How does estimation accuracy improve when additional sensors
(pH, dissolved CO₂, temperature) are incorporated?

**H3 — Fault identification:**
Can EKF residuals not only detect faults but also identify their
root cause through structured residual analysis?

**H4 — Sensor degradation robustness:**
How robust is the EKF against increasing measurement noise (H4a)
and systematic sensor drift over time (H4b)?

---

## Versioning

Each research stage is tagged in git:

| Tag | Description |
|---|---|
| `v1.0-baseline` | Reproduction of Yousefi-Darani et al. (2020) |
| `v2.0-H4a` | Sensor noise robustness experiment |
| `v3.0-H4b` | Sensor drift experiment |
| `v4.0-H2` | Additional sensor inputs |
| `v5.0-H1` | Observability limit analysis |
| `v6.0-H3` | Fault identification |

---

## Requirements

```
python >= 3.10
numpy
scipy
sympy
matplotlib
jupyter
```

Install dependencies:

```bash
pip install numpy scipy sympy matplotlib jupyter
```

---

## Running the notebook

```bash
jupyter notebook ekf_bioreactor.ipynb
```

Run cells sequentially from top to bottom. Each module is
self-contained with explanatory markdown cells.

---

## Why Python instead of MATLAB?

The original paper provides MATLAB code in the appendix. This
implementation uses Python because:

- Free and open source — no license required
- `scipy.integrate.solve_ivp` is equivalent to MATLAB's `ode45`
  (both implement the Dormand-Prince RK45 algorithm)
- `sympy` provides symbolic differentiation equivalent to MATLAB's
  Symbolic Math Toolbox
- Fully reproducible and shareable without software dependencies

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Citation

If you use this code, please cite the original paper:

```bibtex
@Inbook{YousefiDarani2020,
  author    = {Yousefi-Darani, Amirhossein
               and Paquet-Durand, Olivier
               and Hitzmann, Bernd},
  title     = {The Kalman Filter for the Supervision
               of Cultivation Processes},
  booktitle = {Digital Twins},
  year      = {2020},
  publisher = {Springer International Publishing},
  address   = {Cham},
  pages     = {93--124},
  series    = {Advances in Biochemical Engineering/Biotechnology},
  volume    = {177},
  doi       = {10.1007/10_2020_145},
  url       = {https://doi.org/10.1007/10_2020_145}
}
```
