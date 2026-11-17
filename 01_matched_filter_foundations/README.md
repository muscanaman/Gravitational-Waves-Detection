# 01 - Matched filter foundations

Matched filtering for gravitational-wave detection, **derived by hand from first principles** and **verified numerically** on a toy chirp buried in white Gaussian noise, followed by a first look at colored noise and the power spectral density (PSD).

**Author of these notes and code:** Muskan Aman Khan

> **Status: in progress.** Derivations (Pieces 1-5) and the white-noise experiments are complete. Colored noise is done up to generating noise with a chosen PSD and checking it (Steps 7.1-7.2). PSD estimation from a single recording (Welch method) and whitening real data are still to do. See the [roadmap](#6-roadmap).
> **Scope:** the signal is a *toy* linear chirp, not a general-relativity waveform, and no LIGO data is used. The goal is to understand and verify the statistics behind matched-filter searches.

---

## 1. The problem in plain words

A detector records data `x(t) = s(t) + n(t)`: a signal `s(t)` plus noise `n(t)`. If the signal shape is known, what is the best linear filter for deciding whether a signal is present? And what changes when the noise is not white, or when the signal has unknown parameters (phase, arrival time, chirp shape)? These notes answer those questions and check each answer with a simulation.

---

## 2. Mathematics covered (notes + code)

Step names follow the study sequence: Pieces 1-5 are the derivation, Piece 6 the white-noise implementation, Piece 7 colored noise in code.

### Pieces 1-2: setup and the SNR of a linear filter

For a filter output $y=\int x(t)\,h(t)\,dt$ and white Gaussian noise with variance $\sigma^2$, the mean and variance of $y$ are

$$\mathbb{E}[y]=\int s(t)h(t)\,dt \quad(\text{signal present}), \qquad \mathrm{Var}(y)=\sigma^2\!\int h(t)^2\,dt \quad(\text{noise only})$$

The variance uses the white-noise autocorrelation $\mathbb{E}[n(t)n(t')]=\sigma^2\delta(t-t')$ to collapse a double integral into a single one. Hence

$$\mathrm{SNR}=\frac{\mathbb{E}[y]^2}{\mathrm{Var}(y)}=\frac{\left(\int s\,h\,dt\right)^2}{\sigma^2\int h^2\,dt}$$

### Piece 3: the optimal filter (the matched filter)

Treat $s$ and $h$ as vectors. Then $\mathrm{SNR}=\dfrac{\|s\|^2}{\sigma^2}\cos^2\theta$, where $\theta$ is the angle between them (the length of $h$ cancels). By Cauchy-Schwarz,

$$\mathrm{SNR}\le\frac{E_s}{\sigma^2},\quad E_s=\int s^2\,dt,\qquad\text{with equality iff } h(t)\propto s(t)$$

The best filter is a copy of the signal, and for white noise only the signal's *energy* matters, not its shape.

### Piece 4: colored noise, the PSD, and whitening

With noise PSD $S_n(f)$ and Fourier transforms $\tilde s,\tilde H$,

$$\mathrm{SNR}=\frac{\left|\int \tilde s(f)\tilde H^*(f)\,df\right|^2}{\int S_n(f)\,|\tilde H(f)|^2\,df}$$

Substituting $A=\tilde s/\sqrt{S_n}$ and $B=\tilde H\sqrt{S_n}$ turns this into the white-noise problem. Cauchy-Schwarz then gives

$$\tilde H(f)\propto\frac{\tilde s(f)}{S_n(f)},\qquad \rho^2_{\max}=\int\frac{|\tilde s(f)|^2}{S_n(f)}\,df$$

In the gravitational-wave literature this appears as $\rho^2=4\int_0^\infty |\tilde h(f)|^2/S_n(f)\,df$ with a *one-sided* $S_n$; conventions differ by constants, not structure.

### Piece 5: unknown parameters, template banks, trial factor

- **Maximum statistic:** compute the filter output for many parameter values and take the maximum.
- **Arrival time** is handled for free by sliding the template (an inverse FFT in practice).
- **Phase** is removed with two quadrature templates: $|z|=\sqrt{a^2+b^2}$, with $a=(x|s_c)$ and $b=(x|s_q)$.
- **Template bank:** a grid of templates whose spacing is set by the **match** (overlap) between neighbors. A minimal match of 0.97 loses at most about 3% SNR, which is about 9% of the detectable volume (volume scales as $\rho^3$).
- **Trial factor:** the maximum over $N$ independent trials needs a higher threshold, $\rho\approx\sqrt{2\ln(N/p)}$.

### Piece 7 (partial): the PSD and colored noise in code

- One-sided PSD estimate: $S(f_k)=\dfrac{2|X_k|^2}{f_s N}$. For white noise of standard deviation $\sigma$ the expected flat level is $2\sigma^2/f_s$.
- A single recording gives a jagged spectrum; averaging many recordings reveals the true shape.
- Noise with a chosen PSD $S(f)$ is made by filtering white noise with gain $g(f)=\sqrt{S(f)/S_{\text{white}}}$, since amplitude gain $g$ multiplies power by $g^2$.

---

## 3. Experiments and results

All numbers below are produced by the scripts in `code/` (fixed random seeds, so they are reproducible).

**Setup:** sampling rate 1024 Hz; chirp rising from 30 Hz to 250 Hz over 1 s; 4 s of white Gaussian noise with $\sigma$ chosen so the optimal SNR is 8; the chirp is secretly injected at 1.5 s.

| Experiment | Result |
|---|---|
| Signal energy / noise level | $E_s=51.69$, $\sigma=0.8987$, optimal SNR = 8.00 |
| Matched filter (template with the correct phase) | peak $|z|=8.19$ at $t=1.500$ s (true arrival 1.5 s); std of the output away from the peak $=0.987$ (theory: 1) |
| Unknown phase, single template ($\phi_0=0,\ \pi/4,\ \pi/2,\ \pi$) | $z$ at 1.5 s $=8.19,\ 5.85,\ 0.19,\ -7.81$, following $\rho\cos\phi_0$ |
| Unknown phase, two quadratures | $|z|=8.20,\ 8.45,\ 8.44,\ 7.82$: the phase dependence is removed |
| Overlap of the two quadrature templates | $-5\times10^{-13}$ (orthogonal to numerical precision) |
| Noise-only tail $P(|z|>\rho)$ vs $e^{-\rho^2/2}$ ($9.2\times10^5$ samples) | $\rho=1$: 0.604 vs 0.607; $\rho=2$: 0.134 vs 0.135; $\rho=3$: 0.0105 vs 0.0111; $\rho=4$: $3.48\times10^{-4}$ vs $3.35\times10^{-4}$; mean $|z|^2=1.988$ (theory 2) |
| Trial factor (max over 3073 lags, 300 noise-only datasets) | median 3.71, 95th percentile 4.52, largest 5.02; effective independent trials $N_{\text{eff}}\approx680$ (vs 3073 lags) |
| Match vs end-frequency offset ($f_1'$ vs true 250 Hz) | 0.990 (252), 0.945 (255), 0.827 (260), 0.627 (270), 0.434 (290), 0.300 (320) |
| Template-bank search, unknown $f_1$ (truth 259.5 Hz, midway between bank points; 18 templates, 7 Hz spacing) | best template 263 Hz, max $|z|=7.86$ at 1.509 s; the exact template gives 7.93; 1% false-alarm threshold with the bank $=5.29$ |
| White-noise PSD estimate | predicted flat level 0.001953; average of 200 recordings 0.001955; relative scatter 0.97 (one recording) vs 0.068 (average of 200) |
| Colored noise vs target PSD ($S=S_{\text{white}}[1+(50\,\text{Hz}/f)^2]$) | estimated/target ratio between 0.987 and 1.012 in all six frequency bands |

**Reading the results.**

- Matched filtering recovers a chirp that is invisible by eye (SNR 8), at the right time.
- With the wrong phase the output follows $\rho\cos\phi_0$; the two-quadrature statistic removes this.
- The measured noise statistics agree with the theory $e^{-\rho^2/2}$; searching over many lags raises the noise-only maximum (about 3.7 median here), which is why detection thresholds sit well above single-lag values.
- In the bank search, a template 3.5 Hz away still recovers about 98% of the SNR (7.86 vs 7.93 for the exact template in this noise realization), and neighboring wrong templates peak at shifted times (arrival time and chirp rate are partly degenerate).
- Colored noise generated by filtering white noise reproduces the target PSD to within about 1.5% per band.

### Figures (`figures/`)

| File | What it shows |
|---|---|
| `01_toy_chirp.png` | The toy chirp: frequency and amplitude rise with time |
| `02_data_vs_hidden_signal.png` | Data $x=s+n$ (top, signal not visible) and the hidden signal (bottom) |
| `03_matched_filter_output.png` | Matched-filter output vs template start time; a single spike at the true arrival 1.5 s |
| `04_matched_filter_zoom.png` | Zoom of $\pm20$ ms around the spike: a positive peak with oscillating side lobes |
| `05_noise_only_tail_rayleigh.png` | Histogram of noise-only $|z|$ against the Rayleigh curve $r\,e^{-r^2/2}$ |
| `06_trial_factor_max_z.png` | Maximum $|z|$ over all lags for noise-only data, with the injected signal ($|z|\approx8$) far to the right |
| `07_psd_white_noise.png` | A pure 100 Hz tone gives one spike; one white-noise spectrum is jagged; the average of 200 is flat |
| `08_colored_noise_psd.png` | White vs colored noise in time, and the estimated PSD of colored noise on top of the target |

---

## 4. Code (`code/`)

| Script | Step(s) | What it does |
|---|---|---|
| `chirp_step1.py` | 6.1 | Builds and plots the toy chirp |
| `chirp_step2.py` | 6.2-6.9 | Buries the chirp in noise; matched filter by sliding the template; the phase problem; two-quadrature statistic; noise-only Monte Carlo; trial factor; match between templates; template-bank search |
| `psd_step1.py` | 7.1 | Spectrum of a tone, of one white-noise recording, and of the average of many |
| `psd_step2.py` | 7.2 | Generates colored noise with a chosen PSD and checks its estimated PSD against the target |

### How to run

```bash
pip install numpy matplotlib
cd code
python chirp_step2.py     # prints all numbers in the table above (opens several plot windows; close each to continue)
python psd_step2.py
```

The scripts print their results and display the plots; the saved images in `figures/` were exported from these plots.

---

## 5. Notes (`notes/Notes_MatchedFilter_Foundations.pdf`)

Handwritten notes, 31 pages.

- **Pages 1-26** (matched-filter study): the detection problem and Gaussian noise; the SNR definition; mean and variance of the filter output, autocorrelation and the delta function; Cauchy-Schwarz and the matched filter; colored noise, the PSD and whitening; unknown parameters, template banks, the trial factor and the classical search pipeline; the coding steps (building the chirp, burying it in noise, sliding the template); colored noise in code (frequency view, PSD, creating colored noise).
- **Pages 27-31** are the start of notes on a *different* topic: machine-learning noise regression in gravitational-wave detectors (witness channels, sampling and aliasing, the 60 Hz line and its sidebands, least-squares noise regression, and why a linear time-invariant filter cannot remove sidebands). They are kept in this PDF for now and may move to their own study folder.

The code covers more than the notes: the experiments from the phase problem onwards (Steps 6.4-6.9) are in the code and results but not yet written up in the notes.

---

## 6. Roadmap

| Step | Topic | Status |
|---|---|---|
| Pieces 1-5 | Derivation: SNR, optimal filter, colored noise and whitening, unknown parameters, trial factor | done (notes) |
| 6.1-6.9 | White-noise implementation: injection, matched filter, phase, Monte Carlo, trial factor, match, template bank | done (code) |
| 7.1-7.2 | PSD concept; generating colored noise with a chosen PSD | done |
| 7.3 | Estimating the PSD from one long recording (Welch method) and why a window is needed | to do |
| 7.4 | Whitening data and rerunning the matched filter on colored noise | to do |

---

## 7. Limitations

- The signal is a toy chirp with one extra unknown parameter (end frequency); real waveforms depend on masses, spins, sky position and more, and real template banks are far larger.
- The noise is white or synthetic colored noise that is stationary and Gaussian. Real detector noise drifts over time and contains non-Gaussian glitches. Real search thresholds are set empirically from time-shifted background data, not from $e^{-\rho^2/2}$.
- $N_{\text{eff}}$ is a single-parameter fit at the median of the maxima: it slightly under-predicts the tail (about 4.36 predicted vs 4.52 observed at the 95th percentile).
- Single detector only; no network coherence.

---

## 8. References

- S. M. Kay, *Fundamentals of Statistical Signal Processing, Volume II: Detection Theory*, Prentice Hall (1998).
- M. Maggiore, *Gravitational Waves, Volume 1: Theory and Experiments*, Oxford University Press (2008).
- B. J. Owen, "Search templates for gravitational waves from inspiraling binaries: Choice of template spacing", Phys. Rev. D 53, 6749 (1996). *(template spacing and minimal match)*
- B. Allen et al., "FINDCHIRP: an algorithm for detection of gravitational waves from inspiraling compact binaries", Phys. Rev. D 85, 122006 (2012).

*Check the citation details of the last two against the original publications before reusing them elsewhere.*
