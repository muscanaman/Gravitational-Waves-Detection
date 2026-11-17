# Gravitational-Waves-Detection

Hand-derived and numerically verified methods for **gravitational-wave (GW) data analysis**: finding signals buried in instrumental noise. Each study works through one method or paper in three stages: **derive the mathematics by hand, check it with a small simulation, and write down the gaps and open questions.**

**Author:** Muskan Aman Khan

> This is a research-preparation and learning repository, written while building a foundation in GW data analysis. Everything here is a study of known methods on toy or public data. Nothing in it claims a new result, and open questions are labeled as questions.

---

## Why this repository

Gravitational-wave detectors record data that is dominated by noise. Detecting a signal means understanding that noise: its spectrum, its non-Gaussian "glitches", and what can be removed without destroying the signal. The studies here follow that thread:

1. **Matched filtering**, the classical optimal detector, and the noise statistics behind it (the baseline everything else is compared against).
2. **Denoising**, removing noise from the data without erasing sharp signal features.
3. **Glitch classification**, identifying non-Gaussian noise transients with a neural network and explaining its decisions.

The goal of each folder is to understand *why* a method works, in the mathematics, and to verify that understanding numerically.

---

## Studies

| Folder | Study | Based on | Status |
|---|---|---|---|
| [`01_matched_filter_foundations`](01_matched_filter_foundations/) | Matched filtering from first principles: optimal SNR, colored noise and whitening, unknown parameters, template banks, false-alarm statistics; toy-chirp experiments; PSD and colored noise | S. M. Kay, *Detection Theory*; M. Maggiore, *Gravitational Waves Vol. 1* | In progress: derivations and white-noise experiments done; PSD estimation and whitening in code to do |
| [`02_tv_denoising_rROF_cWB`](02_tv_denoising_rROF_cWB/) | Total-variation (rROF) denoising and its use in the cWB search pipeline: classical smoothing and why it blurs edges, then total variation | P. J. Barneo et al., arXiv:2202.10465 (2022) | In progress: first two steps done; total variation paused |
| [`03_glitch_classification_scoreCAM_UMAP`](03_glitch_classification_scoreCAM_UMAP/) | Glitch classification with a CNN, explained with Score-CAM and UMAP | N. Koyama et al., Mach. Learn.: Sci. Technol. 5, 035028 (2024) | In progress: notes under way, code not started |

Planned, not started: a neural-network detector compared with the matched filter on identical data, and noise regression with witness channels.

**Suggested reading order:** 01 (the noise model and the optimal baseline) → 02 (removing noise) → 03 (non-Gaussian noise).

---

## How each folder is organized

```
NN_study_name/
├── README.md     what the study is, the math covered, results, limitations, next steps
├── code/         Python scripts for the experiments
├── notes/        handwritten or written study notes (PDF)
├── figures/      plots produced by the code
└── paper/        paper PDF, where the license allows (study 03 only)
```

Each study's README states exactly what has been done, what the numbers are, and what is still missing.

---

## Getting started

```bash
git clone https://github.com/muscanaman/Gravitational-Waves-Detection.git
cd Gravitational-Waves-Detection
pip install -r requirements.txt
```

Studies 01 and 02 need only NumPy and Matplotlib. For example:

```bash
cd 01_matched_filter_foundations/code
python chirp_step2.py        # white-noise matched filtering experiments
```

Study 03 will have its own requirements when its code is added. Each study's README says how to run its scripts.

---

## Scope and honesty notes

- The experiments use synthetic signals and noise (a toy chirp, a noisy box signal). They illustrate the mathematics; they are not reproductions of published results on LIGO data.
- Real detector noise is non-stationary and non-Gaussian; the simple statistical results derived here (for example, false-alarm formulas) hold for idealized noise, and each study states its assumptions.
- Paper PDFs are included only where the license allows it (study 03's paper is open access under CC BY 4.0). Other papers are cited by reference only.

---

## Author

**Muskan Aman Khan**: electrical engineering graduate (Namal University), focused on the mathematics of machine learning and its application to gravitational-wave data analysis.
