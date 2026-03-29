# Computational Simulation of Drug Diffusion Through Biological Tissue
### Modeling Fick's Second Law, Layered Tissue, and Molecular Binding Kinetics in Python

**Author:** Abhir Mahajan | **School:** Lynbrook High School | **Year:** 2028  
**Target Competition:** Synopsys Championship (Santa Clara County ISEF Affiliate)

---

## Overview

This project computationally simulates how drugs diffuse through biological tissue using Fick's second law of diffusion. The goal is to model the physical and biochemical factors that determine how much of a drug reaches its target site — a question with direct relevance to drug dosing, therapeutic effectiveness, and pharmaceutical design.

Rather than observing drug behavior clinically from the outside, this simulation models diffusion at the molecular level: how concentration gradients, tissue geometry, time, and binding interactions shape the drug's journey through tissue.

---

## Motivation

In a clinical setting, a physician sees outcomes — did the drug work or not? What they cannot easily see is *why*: how the drug moved through tissue, where it was lost to binding, and at what depth the concentration became therapeutically ineffective.

Computational simulation makes these invisible processes visible. This project is a step toward building models that can predict therapeutic effectiveness from first principles, before a drug ever reaches a patient.

---

## Experiments

| # | Condition Tested | Key Finding |
|---|-----------------|-------------|
| 1 | Baseline diffusion | Established reference concentration profile across tissue depth |
| 2 | Diffusion across multiple time points | Concentration gradients flatten and penetrate deeper over time |
| 3 | Varying diffusion coefficient (D) | Higher D accelerates penetration; small changes in D produce large depth effects |
| 4 | Layered tissue (heterogeneous medium) | Tissue boundaries create discontinuities in concentration profiles, altering drug delivery to deeper layers |
| 5 | Diffusion with molecular binding | Binding substantially reduces free drug concentration and limits penetration depth across all tissue layers |

---

## Key Result: Experiment 5

The graph below compares normal diffusion against diffusion with a binding term:

![Diffusion with Binding](figures/experiment5_binding.png)

**Observations:**
- Free drug concentration decreases more rapidly in the presence of binding at all tissue depths
- Overall drug magnitude is reduced across the entire tissue profile
- Effective penetration depth is smaller when binding effects are present
- The binding effect compounds over time

**Interpretation:**  
Biochemical interactions within tissue — such as binding to proteins or cellular uptake — substantially reduce the amount of free drug available to diffuse. This highlights the importance of accounting for drug-tissue interactions when predicting therapeutic effectiveness.

---

## Current Limitation: Binding Model (Important)

**This section is included in the interest of scientific transparency.**

The binding term in Experiment 5 is currently implemented as a first-order decay multiplied onto the analytical diffusion solution:

```python
C_bind = np.exp(-(x**2) / (4 * D * t)) * np.exp(-k_bind * t)
```

This models uniform exponential decay over time — it is a valid approximation of drug degradation or nonspecific uptake, but it is **not** a physically accurate model of reversible molecular binding. True binding kinetics are:
- Concentration-dependent
- Saturable (binding sites get used up)
- Reversible (drug can unbind)

### Planned Implementation

The next version will replace this with a coupled PDE system:

```
∂C/∂t       = D(∂²C/∂x²) - k_on·C·B + k_off·C_bound
∂C_bound/∂t = k_on·C·B - k_off·C_bound
∂B/∂t       = -k_on·C·B + k_off·C_bound
```

Where `C` is free drug concentration, `C_bound` is bound drug, `B` is available binding site concentration, and `k_on`/`k_off` are binding/unbinding rate constants drawn from pharmacokinetic literature.

This implementation is currently in progress.

---

## Repository Structure

```
/
├── README.md
├── experiment1_baseline.py
├── experiment2_time.py
├── experiment3_diffusion_coefficient.py
├── experiment4_layered_tissue.py
├── experiment5_binding.py
├── figures/
│   ├── experiment1_baseline.png
│   ├── experiment2_time.png
│   ├── experiment3_D_variation.png
│   ├── experiment4_layered.png
│   └── experiment5_binding.png
└── analysis/
    └── final_analysis.md
```

---

## How to Run

**Requirements:**
```
Python 3.8+
numpy
matplotlib
```

**Install dependencies:**
```bash
pip install numpy matplotlib
```

**Run any experiment:**
```bash
python experiment1_baseline.py
python experiment5_binding.py
# etc.
```

Each script generates and displays its corresponding figure.

---

## Next Steps

- [ ] Implement reversible binding kinetics via coupled PDE system (Experiment 5 revision)
- [ ] Incorporate biologically validated parameters (D, k_on, k_off) from pharmacokinetic literature for a specific drug class
- [ ] Extend layered tissue model (Experiment 4) to use tissue-specific diffusion coefficients
- [ ] Combine layered tissue + binding kinetics into a unified model
- [ ] Explore numerical stability tradeoffs between explicit Euler and Crank-Nicolson methods

---

## Contact

Abhir Mahajan
Lynbrook High School
abhirymahajan@gmail.com
