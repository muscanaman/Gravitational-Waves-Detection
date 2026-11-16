import numpy as np
import matplotlib.pyplot as plt

fs = 1024                    # samples per second
T = 4.0                      # seconds of data
N = int(T * fs)              # number of samples
t = np.arange(N) / fs
rng = np.random.default_rng(0)

def psd_one_sided(x, fs):
    """Raw spectrum estimate: how much power sits at each frequency."""
    X = np.fft.rfft(x)                        # FFT for real data (positive frequencies)
    f = np.fft.rfftfreq(len(x), 1 / fs)       # the frequency belonging to each bin
    P = 2 * np.abs(X)**2 / (fs * len(x))      # the formula above
    return f, P

# --- Demo 1: a pure 100 Hz tone ---
sine = np.sin(2 * np.pi * 100 * t)
f, P_sine = psd_one_sided(sine, fs)
print("frequency spacing:", f[1] - f[0], " number of bins:", len(f))
print("sine peak at f =", f[np.argmax(P_sine)], "Hz")

# --- Demo 2: white noise, one recording vs the average of many ---
sigma = 1.0
f, P1 = psd_one_sided(sigma * rng.standard_normal(N), fs)
M = 200
Pavg = np.mean([psd_one_sided(sigma * rng.standard_normal(N), fs)[1] for _ in range(M)], axis=0)

print("predicted flat level 2*sigma^2/fs =", 2 * sigma**2 / fs)
print("single recording: mean =", P1[1:-1].mean(), " spread/mean =", P1[1:-1].std() / P1[1:-1].mean())
print("average of", M, ":   mean =", Pavg[1:-1].mean(), " spread/mean =", Pavg[1:-1].std() / Pavg[1:-1].mean())

fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(f, P_sine); ax[0].set_title("100 Hz sine: one spike")
ax[1].plot(f, P1, lw=0.5); ax[1].set_title("White noise: ONE recording (jagged)")
ax[2].plot(f, Pavg, lw=0.8); ax[2].axhline(2 * sigma**2 / fs, color="r", ls="--", label="theory")
ax[2].set_title(f"White noise: average of {M} (flat)"); ax[2].legend()
for a in ax: a.set_xlabel("frequency (Hz)")
ax[0].set_ylabel("power per Hz")
plt.tight_layout(); plt.show()