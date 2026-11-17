# 03 - Glitch Classification with CNN, Score-CAM and UMAP

Study of: Koyama N. et al., "Enhancing the rationale of convolutional neural networks
for glitch classification in gravitational wave detectors: a visual explanation",
Mach. Learn.: Sci. Technol. 5, 035028 (2024). https://doi.org/10.1088/2632-2153/ad6391
(Open access, CC BY 4.0. Paper PDF stored in `paper/` with attribution.)

## Why this paper
Glitches (non-stationary, non-Gaussian detector noise) cause false alarms and bias
parameter estimation. This paper classifies 22 Gravity Spy glitch classes with a CNN
and explains the decisions with Score-CAM saliency maps and UMAP feature-space plots.

## Goals of this folder
1. Written notes covering the paper's mathematics, results, gaps and open questions.
2. A from-scratch reproduction of the pipeline (CNN, Score-CAM, UMAP) in PyTorch.
3. Extensions that address gaps in the paper (ablation, imbalance, calibration,
   Score-CAM faithfulness, UMAP stability).

## Structure

    03_glitch_classification_scoreCAM_UMAP/
        code/       Python scripts
        figures/    generated plots
        notes/      study notes (markdown)
        paper/      paper PDF
        data/       Gravity Spy data (not committed, see .gitignore)

## Status
- [x] Scaffold and paper PDF
- [ ] Notes: written through the UMAP cost function (Eq. 7); remaining: Eq. 8, results and evaluation, gaps
- [ ] Data download and inspection
- [ ] CNN training and evaluation
- [ ] Score-CAM
- [ ] UMAP
- [ ] Extensions
