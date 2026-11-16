import numpy as np
import matplotlib.pyplot as plt

# ---- Step 6.1: the chirp (same as before) ----
fs = 1024
D = 1.0
f0, f1 = 30.0, 250.0
tau = np.arange(0, D, 1/fs)
k = (f1 - f0) / D
phase = 2*np.pi*(f0*tau + 0.5*k*tau**2)
s = (tau / D) * np.hanning(len(tau)) * np.sin(phase)

# ---- Step 6.2: hide it in noise ----
rng = np.random.default_rng(0)       # fixed seed -> reproducible
T = 4.0                              # total data length in seconds
N = int(T * fs)                      # number of samples = 4096
t = np.arange(N) / fs                # time axis for the long stretch

t_inject = 1.5                       # secret arrival time
i0 = int(t_inject * fs)              # same moment, as an array index
signal = np.zeros(N)
signal[i0:i0 + len(s)] = s           # chirp placed inside a stretch of zeros

E_s = np.sum(s**2)                   # signal energy
rho_target = 8.0
sigma = np.sqrt(E_s) / rho_target    # noise level giving SNR = 8
noise = sigma * rng.standard_normal(N)

x = signal + noise                   # what the "detector" records

print("signal energy:", E_s)
print("sigma:", sigma)
print("optimal SNR (check):", np.sqrt(E_s) / sigma)

fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
ax[0].plot(t, x, lw=0.6);      ax[0].set_title("Data x = s + n  (what we actually see)")
ax[1].plot(t, signal, lw=0.8); ax[1].set_title("Hidden signal s (truth, we normally don't know this)")
ax[1].set_xlabel("time (s)")
plt.tight_layout(); plt.show()

# ---- Step 6.3: matched filter by sliding the template ----
# np.correlate(x, s, 'valid')[j] = sum_i x[j+i] * s[i]
z = np.correlate(x, s, mode='valid') / (sigma * np.sqrt(E_s))

t_z = np.arange(len(z)) / fs        # z[j] <-> template STARTING at time j/fs

j_peak = np.argmax(np.abs(z))
print("peak |z| =", np.abs(z[j_peak]))
print("peak time =", t_z[j_peak], "s   (true arrival = 1.5 s)")
print("std of z away from the peak:", np.std(np.delete(z, np.arange(j_peak-50, j_peak+50))))

plt.figure(figsize=(10, 4))
plt.plot(t_z, z, lw=0.7)
plt.axvline(t_inject, color='r', ls='--', label='true arrival (1.5 s)')
plt.xlabel("template start time (s)"); plt.ylabel("matched-filter output z")
plt.title("Matched filter output (units of noise std)")
plt.legend(); plt.tight_layout(); plt.show()

mask = (t_z > 1.48) & (t_z < 1.52)
plt.figure(figsize=(8, 3))
plt.plot(t_z[mask], z[mask], 'o-', ms=3)
plt.axvline(t_inject, color='r', ls='--')
plt.xlabel("template start time (s)"); plt.ylabel("z")
plt.title("Zoom: ±20 ms around true arrival")
plt.tight_layout(); plt.show()

print("z at exactly 1.5 s:", z[int(t_inject * fs)])

# ---- Step 6.4: what happens if the signal's phase is not what the template assumes? ----
amp = (tau / D) * np.hanning(len(tau))

for phi0 in [0, np.pi/4, np.pi/2, np.pi]:
    s_phi = amp * np.sin(phase + phi0)           # signal with an extra phase offset
    sig = np.zeros(N); sig[i0:i0 + len(s)] = s_phi
    x_phi = sig + noise                          # SAME noise as before
    z_phi = np.correlate(x_phi, s, mode='valid') / (sigma * np.sqrt(E_s))   # template still phase 0
    print(f"phi0 = {phi0:4.2f}   z at 1.5 s = {z_phi[i0]:6.2f}   "
          f"predicted 8*cos(phi0) = {8*np.cos(phi0):6.2f}   max|z| = {np.max(np.abs(z_phi)):.2f}")
    
    # ---- Step 6.5: two-quadrature filter ----
s_c = amp * np.sin(phase)            # phase-0 template (same as s)
s_q = amp * np.cos(phase)            # same chirp, shifted by 90 degrees
E_c, E_q = np.sum(s_c**2), np.sum(s_q**2)
print("overlap (s_c . s_q)/sqrt(E_c E_q) =", np.sum(s_c * s_q) / np.sqrt(E_c * E_q))

def mf_mag(xdata):
    a = np.correlate(xdata, s_c, 'valid') / (sigma * np.sqrt(E_c))
    b = np.correlate(xdata, s_q, 'valid') / (sigma * np.sqrt(E_q))
    return np.sqrt(a**2 + b**2)

for phi0 in [0, np.pi/4, np.pi/2, np.pi]:
    sig = np.zeros(N); sig[i0:i0 + len(s)] = amp * np.sin(phase + phi0)
    m = mf_mag(sig + noise)                      # same noise again
    jm = np.argmax(m)
    print(f"phi0 = {phi0:4.2f}   |z| at 1.5 s = {m[i0]:5.2f}   "
          f"peak |z| = {m[jm]:5.2f} at t = {t_z[jm]:.3f} s")
    
    # ---- Step 6.6: noise-only Monte Carlo ----
rng2 = np.random.default_rng(1)
M = 300                                   # number of noise-only datasets
mags = []
for _ in range(M):
    n_only = sigma * rng2.standard_normal(N)
    mags.append(mf_mag(n_only))
mags = np.concatenate(mags)               # |z| at every lag of every trial

print("number of |z| samples:", mags.size)
print("mean of |z|^2 (theory 2):", np.mean(mags**2))
for rho in [1, 2, 3, 4]:
    emp = np.mean(mags > rho)
    print(f"rho={rho}: empirical P(|z|>rho) = {emp:.2e}   theory e^(-rho^2/2) = {np.exp(-rho**2/2):.2e}")

plt.figure(figsize=(7, 4))
plt.hist(mags, bins=60, density=True, alpha=0.6, label="simulation")
r = np.linspace(0, 5, 300)
plt.plot(r, r * np.exp(-r**2 / 2), 'r', label="theory: r e^{-r²/2}")
plt.xlabel("|z| (noise only)"); plt.ylabel("density"); plt.legend()
plt.tight_layout(); plt.show()

# ---- Step 6.7: trial factor ----
maxes = mags.reshape(M, -1).max(axis=1)    # one maximum per noise-only dataset
L = mags.size // M                         # lags per dataset (should be 3073)
print("lags per dataset L =", L)
print("median of max|z| (noise only):", np.median(maxes))
print("95th percentile of max|z|:", np.percentile(maxes, 95))
print("largest max in", M, "datasets:", maxes.max())

# prediction if all L lags were independent: median solves 1-(1-p)^L = 0.5
p_single = 1 - 0.5**(1 / L)
print("naive median prediction (L independent lags):", np.sqrt(-2 * np.log(p_single)))

# effective number of independent trials implied by the simulated median
med = np.median(maxes)
N_eff = np.log(0.5) / np.log(1 - np.exp(-med**2 / 2))
print("effective independent trials N_eff ~", N_eff, "  (vs L =", L, ")")

plt.figure(figsize=(7, 4))
plt.hist(maxes, bins=25, alpha=0.7)
plt.axvline(8, color='r', ls='--', label='our injected signal (|z| ~ 8)')
plt.xlabel("max |z| over all lags, noise only"); plt.ylabel("count"); plt.legend()
plt.tight_layout(); plt.show()

# ---- Step 6.8: match between templates with different f1 ----
def make_template(f1_, phase0=0.0):
    k_ = (f1_ - f0) / D
    ph = 2*np.pi*(f0*tau + 0.5*k_*tau**2)
    return amp * np.sin(ph + phase0)

A_c = make_template(250.0, 0.0)             # reference chirp (true f1 = 250)
A_q = make_template(250.0, np.pi/2)         # its 90-degree-shifted version

for f1_ in [250, 252, 255, 260, 270, 290, 320]:
    B = make_template(f1_)                  # trial template with a different f1
    c = np.correlate(B, A_c, mode='full')   # all relative time shifts
    q = np.correlate(B, A_q, mode='full')
    m = np.max(np.sqrt(c**2 + q**2)) / (np.linalg.norm(A_c) * np.linalg.norm(B))
    print(f"f1' = {f1_:4d} Hz   match = {m:.3f}   SNR kept = {100*m:.1f}%")
    
    # ---- Step 6.9: search with unknown f1 using a template bank ----
def mf_mag_f1(xdata, f1_):
    sc = make_template(f1_, 0.0)
    sq = make_template(f1_, np.pi/2)
    a = np.correlate(xdata, sc, 'valid') / (sigma * np.linalg.norm(sc))
    b = np.correlate(xdata, sq, 'valid') / (sigma * np.linalg.norm(sq))
    return np.sqrt(a**2 + b**2)

f1_true, phi_true = 259.5, 1.0                  # secret truth, worst-case between grid points
sig_u = np.zeros(N); sig_u[i0:i0 + len(s)] = make_template(f1_true, phi_true)
x_u = sig_u + noise
print("true optimal SNR:", np.linalg.norm(sig_u) / sigma)

bank = np.arange(200.0, 321.0, 7.0)              # f1 grid, spacing 7 Hz
best = (0.0, None, None)
for f1_ in bank:
    m = mf_mag_f1(x_u, f1_)
    j = np.argmax(m)
    print(f"template f1 = {f1_:6.1f}   max|z| = {m[j]:5.2f} at t = {t_z[j]:.3f} s")
    if m[j] > best[0]:
        best = (m[j], f1_, t_z[j])
print("BEST (|z|, f1, time):", best)
print("with the exact template (not available in practice):", mf_mag_f1(x_u, f1_true).max())

# 1% false-alarm threshold when we also trial the whole bank
N_eff, n_bank = 680, len(bank)
p_single = 1 - 0.99**(1 / (N_eff * n_bank))
print("1% false-alarm threshold with bank:", np.sqrt(-2 * np.log(p_single)))