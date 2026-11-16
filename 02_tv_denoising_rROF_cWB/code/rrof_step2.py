## Step 2.1: create the clean signal and the noisy data
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
N = 400
t = np.arange(N)
u_true = np.where((t > 150) & (t < 250), 1.0, 0.0)
sigma = 0.3
f = u_true + sigma * rng.standard_normal(N)

print("length of f:", len(f))
print("number of samples inside the box:", int(u_true.sum()))

plt.plot(t, f, lw=0.6, color="gray", label="data f = u + n")
plt.plot(t, u_true, "k", lw=1.5, label="true u")
plt.legend(); plt.show()

## Step 2.2: build the “neighbor jump” tool $D$ and see what it measures.

D = np.zeros((N - 1, N))
for i in range(N - 1):
    D[i, i], D[i, i + 1] = -1, 1

print("shape of D:", D.shape)

jumps_true = D @ u_true
print("positions where the true box jumps:", np.nonzero(jumps_true)[0], " values:", jumps_true[np.nonzero(jumps_true)[0]])

print("R(u_true) = ||D u_true||^2 =", np.sum(jumps_true**2))
print("R(f)      = ||D f||^2      =", np.sum((D @ f)**2))

DtD = D.T @ D

def smooth(f, lam):
    """Solve (I + (2/lam) D^T D) u = f."""
    return np.linalg.solve(np.eye(N) + (2 / lam) * DtD, f)

print("RMSE of raw data f:", np.sqrt(np.mean((f - u_true)**2)))
for lam in [0.01, 0.1, 0.5, 1, 5, 20, 100, 1000]:
    u = smooth(f, lam)
    print(f"lambda = {lam:7.2f}   RMSE = {np.sqrt(np.mean((u - u_true)**2)):.3f}")
    import os

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
ax[0].plot(t, f, lw=0.5, color="gray", label="data f")
ax[0].plot(t, u_true, "k", lw=1.5, label="true u")
ax[0].plot(t, smooth(f, 0.1), "r", lw=1.5, label="lambda = 0.1 (best RMSE)")
ax[0].plot(t, smooth(f, 5), "b", lw=1.2, label="lambda = 5")
ax[0].legend(); ax[0].set_title("Smoothness denoising")

nu = np.linspace(0, 0.5, 400)                       # frequency in cycles per sample
for lam in [0.1, 1, 10]:
    ax[1].plot(nu, 1 / (1 + (2 / lam) * 4 * np.sin(np.pi * nu)**2), label=f"lambda = {lam}")
ax[1].set_xlabel("frequency (cycles/sample)"); ax[1].set_ylabel("gain H")
ax[1].legend(); ax[1].set_title("The filter is a low-pass")
plt.tight_layout()

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures", "step2_smoothing.png")
plt.savefig(out, dpi=150)
plt.show()