import numpy as np
import matplotlib.pyplot as plt

fs = 1024 # Sampling frequency
D = 1.0 # Duration of the signal in seconds
f0,f1 = 30.0, 250.0  # Start and end frequencies of the chirp signal


tau = np.arange(0, D, 1/fs) # Time vector
k = (f1 - f0)/D # Chirp rate. How fast the frequency changes with time
phase = 2*np.pi*(f0*tau + 0.5*k*tau**2) # Phase of the chirp signal
amplitude = (tau/D) * np.hanning(len(tau)) # Amplitude envelope of the chirp signal
s = amplitude * np.sin(phase) # Chirp signal

plt.plot(tau, s)
plt.xlabel('Time [s]')
plt.ylabel('Amplitude, s(t)')
plt.title('Toy chirp: frequency and amplitude rise')
plt.show()