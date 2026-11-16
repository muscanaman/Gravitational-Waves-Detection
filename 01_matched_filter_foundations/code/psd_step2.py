import numpy as np
import matplotlib.pyplot as plt

fs = 1024
T = 4.0
N = int(T * fs)
t = np.arange(N) / fs
rng = np.random.default_rng(0)

def psd_one_sided(x, fs):
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), 1 / fs)
    return f, 2 * np.abs(X)**2 / (fs * len(x))

# --- the target noise "fingerprint" ---
f = np.fft.rfftfreq(N, 1 / fs)
S_white = 2 / fs                                  # PSD of white noise with sigma = 1
f_knee = 50.0
shape = 1 + (f_knee / np.maximum(f, 1.0))**2      # S(f) / S_white
S_target = S_white * shape
gain = np.sqrt(shape)                             # amplitude gain = sqrt(power gain)
gain[0] = 0.0                                     # remove the constant (0 Hz) offset

def colored_noise(rng):
    w = rng.standard_normal(N)                    # step 1: white noise
    W = np.fft.rfft(w)                            # step 2: to frequency view
    return np.fft.irfft(W * gain, n=N)            # step 3: reshape, go back to time

# --- check: does the PSD of our noise match the target? ---
M = 200
Pavg = np.mean([psd_one_sided(colored_noise(rng), fs)[1] for _ in range(M)], axis=0)

print("--- ratio (estimated PSD / target PSD), averaged in frequency bands ---")
for lo, hi in [(1, 5), (5, 20), (20, 60), (60, 150), (150, 300), (300, 500)]:
    m = (f >= lo) & (f < hi)
    print(f"{lo:4d}-{hi:4d} Hz: ratio = {Pavg[m].mean() / S_target[m].mean():.3f}")

x_col = colored_noise(rng)
x_white = rng.standard_normal(N)
print("std white:", x_white.std(), "  std colored:", x_col.std())

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
ax[0].plot(t[:1024], x_white[:1024], lw=0.6, label="white")
ax[0].plot(t[:1024], x_col[:1024], lw=0.8, label="colored")
ax[0].set_xlabel("time (s)"); ax[0].set_title("First second of each"); ax[0].legend()
ax[1].loglog(f[1:], Pavg[1:], label="estimated (average of 200)")
ax[1].loglog(f[1:], S_target[1:], "r--", label="target")
ax[1].set_xlabel("frequency (Hz)"); ax[1].set_ylabel("PSD"); ax[1].legend()
plt.tight_layout(); plt.show()