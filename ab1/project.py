import numpy as np
import matplotlib.pyplot as plt

def DFT(x):
    N = len(x)
    X = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            X[k] += x[n] * np.exp(-2j * np.pi * k * n / N)
    return X

def FFT_radix2(x):
    N = len(x)
    if N <= 1:
        return x
    elif N % 2 > 0:
        raise ValueError("N deve ser uma potência de 2")
    even = FFT_radix2(x[0::2])
    odd = FFT_radix2(x[1::2])
    factor = np.exp(-2j * np.pi * np.arange(N) / N)
    return np.concatenate([even + factor[:N//2] * odd,
                           even + factor[N//2:] * odd])

def FFT_radix3(x):
    N = len(x)
    if N == 1:
        return x
    elif N % 3 > 0:
        raise ValueError("N deve ser uma potência de 3")
    x0 = FFT_radix3(x[0::3])
    x1 = FFT_radix3(x[1::3])
    x2 = FFT_radix3(x[2::3])
    factor = np.exp(-2j * np.pi * np.arange(N) / N)
    return np.concatenate([x0 + factor[:N//3] * x1 + factor[:N//3]**2 * x2,
                           x0 + factor[N//3:2*N//3] * x1 + factor[N//3:2*N//3]**2 * x2,
                           x0 + factor[2*N//3:] * x1 + factor[2*N//3:]**2 * x2])

# Escolher um N compatível com radix-2 e radix-3
N_radix2 = 8  # Potência de 2
N_radix3 = 9  # Potência de 3

# Geração do sinal de teste
fs = 100  # Frequência de amostragem
f1, f2 = 10, 20  # Frequências do sinal
n_radix2 = np.arange(N_radix2)
n_radix3 = np.arange(N_radix3)
signal_radix2 = np.sin(2 * np.pi * f1 * n_radix2 / fs) + np.sin(2 * np.pi * f2 * n_radix2 / fs)
signal_radix3 = np.sin(2 * np.pi * f1 * n_radix3 / fs) + np.sin(2 * np.pi * f2 * n_radix3 / fs)

# Aplicação das transformadas
X_dft2 = DFT(signal_radix2)  # DFT para radix-2
X_dft3 = DFT(signal_radix3)  # DFT para radix-3
X_fft2 = FFT_radix2(signal_radix2)
X_fft3 = FFT_radix3(signal_radix3)
X_fft_np2 = np.fft.fft(signal_radix2)  # FFT NumPy para N=8
X_fft_np3 = np.fft.fft(signal_radix3)  # FFT NumPy para N=9

# Plot dos resultados
plt.figure(figsize=(18, 12))
plt.subplot(3, 2, 1)
plt.stem(np.abs(X_dft2))
plt.title("DFT (N=8)")
plt.subplot(3, 2, 2)
plt.stem(np.abs(X_fft2))
plt.title("Radix-2 FFT (N=8)")
plt.subplot(3, 2, 3)
plt.stem(np.abs(X_dft3))
plt.title("DFT (N=9)")
plt.subplot(3, 2, 4)
plt.stem(np.abs(X_fft3))
plt.title("Radix-3 FFT (N=9)")
plt.subplot(3, 2, 5)
plt.stem(np.abs(X_fft_np2))
plt.title("NumPy FFT (N=8)")
plt.subplot(3, 2, 6)
plt.stem(np.abs(X_fft_np3))
plt.title("NumPy FFT (N=9)")
plt.show()
