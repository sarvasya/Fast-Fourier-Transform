# 📡 FFT Demo — C++ Digital Signal Processing

> Cooley-Tukey Radix-2 FFT implementation in pure C++17 — no external dependencies.

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/)
[![License: Public Domain](https://img.shields.io/badge/license-Public%20Domain-green.svg)]()
[![Build: GCC/Clang/MSVC](https://img.shields.io/badge/build-GCC%20%7C%20Clang%20%7C%20MSVC-lightgrey.svg)]()

---

## 📖 Overview

This program demonstrates the **Fast Fourier Transform (FFT)** and its core applications in **Digital Signal Processing (DSP)**:

- Synthesises a composite multi-tone signal in the time domain
- Transforms it to the frequency domain using the **Cooley-Tukey Radix-2 DIT FFT**
- Analyses the **magnitude** and **phase** spectrum
- Reconstructs the original signal via the **Inverse FFT (IFFT)**
- Applies a **frequency-domain low-pass filter**

All in ~250 lines of standard C++ using only `<complex>`, `<cmath>`, and the STL.

---

## 🧮 The Math

### Discrete Fourier Transform (DFT)

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j \frac{2\pi k n}{N}}$$

A naive DFT runs in **O(N²)**. The FFT reduces this to **O(N log₂ N)** using the Cooley-Tukey divide-and-conquer approach.

### Butterfly Operation

Each FFT stage combines pairs of values using:

```
top    = E[k] + W^k · O[k]
bottom = E[k] − W^k · O[k]
```

where `W = e^{-j 2π/N}` is the **twiddle factor** (principal root of unity).

---

## ⚡ Quick Start

### Prerequisites

- C++17 compiler: `g++ 7+`, `clang++ 5+`, or `MSVC 2017+`
- No external libraries required

### Build & Run

**Linux / macOS**
```bash
g++ -std=c++17 -O2 -o fft_demo fft_demo.cpp
./fft_demo
```

**Windows (MSVC)**
```bat
cl /std:c++17 /O2 fft_demo.cpp /Fe:fft_demo.exe
fft_demo.exe
```

**Windows (MinGW / MSYS2)**
```bash
g++ -std=c++17 -O2 -o fft_demo.exe fft_demo.cpp
./fft_demo.exe
```

---

## 🔬 What the Demo Does

| Step | Description |
|------|-------------|
| **1. Signal synthesis** | Generates a composite sine wave: 50 Hz (amp 1.0) + 120 Hz (amp 0.5) + 200 Hz (amp 0.75) |
| **2. Forward FFT** | Transforms the signal to the frequency domain in-place using Cooley-Tukey |
| **3. Spectrum analysis** | Computes magnitude & phase; identifies peaks; prints an ASCII bar chart |
| **4. Inverse FFT** | Reconstructs the original signal; verifies round-trip error ≈ machine epsilon |
| **5. Low-pass filter** | Zeros all bins above 100 Hz; only the 50 Hz tone survives |

---

## 🖥️ Sample Output

```
╔══════════════════════════════════════════════════╗
║   FFT Demo — Digital Signal Processing in C++   ║
╚══════════════════════════════════════════════════╝

Signal composition:
  Tone A: 50 Hz, amplitude 1
  Tone B: 120 Hz, amplitude 0.5
  Tone C: 200 Hz, amplitude 0.75

── Frequency spectrum (peaks > 0.1 amplitude) ──
  f =   50.00 Hz   |amplitude| = 1.0000   phase = -1.5708 rad
  f =  200.00 Hz   |amplitude| = 0.7500   phase = -1.5708 rad
  f =  120.00 Hz   |amplitude| = 0.5000   phase = -1.5708 rad

── Magnitude spectrum (one-sided) ──
       50 Hz │████████████████████████████████████ 1.0
      120 Hz │██████████████████ 0.5
      200 Hz │███████████████████████████ 0.75

── Inverse FFT (IFFT) reconstruction ──
  Max round-trip error: 2.78e-14  (near machine epsilon)

── Low-pass filter (cutoff = 100 Hz) ──
  Dominant frequency after filtering: 50.0 Hz  (expected 50 Hz)

✓ Demo complete.
```

---

## 🗂️ Code Structure

```
fft_demo.cpp
│
├── isPowerOfTwo(n)        — validates FFT size constraint
├── zeroPad(x)             — pads signal to next power of 2
├── fft(x, inverse)        — in-place Cooley-Tukey FFT / IFFT
├── magnitudeSpectrum(X)   — returns |X[k]| for one-sided spectrum
├── phaseSpectrum(X)       — returns arg(X[k]) in radians
├── binToHz(k, N, fs)      — converts bin index → frequency in Hz
├── printBar(...)          — ASCII bar chart renderer
└── main()                 — orchestrates all five demo steps
```

---

## 📐 Key DSP Concepts Demonstrated

| Concept | Where in code |
|---------|--------------|
| Bit-reversal permutation | `fft()` — reorders samples before butterfly stages |
| Twiddle factors | `fft()` — roots of unity `W = e^{-j2π/N}` |
| One-sided spectrum | `magnitudeSpectrum()` — bins 0 to N/2 |
| Amplitude normalisation | `main()` — factor `2/N` for single-sided peak amplitude |
| Rectangular window LPF | `main()` step 5 — zero out bins above cutoff |
| Round-trip fidelity | `main()` step 4 — IFFT ∘ FFT ≈ identity |

---

## 🚀 Extensions

Ideas for taking this further:

- **Windowing** — apply a Hann or Hamming window before the FFT to reduce spectral leakage
- **STFT** — slide a windowed FFT across a longer signal to produce a spectrogram
- **Convolution** — multiply two FFTs element-wise and IFFT for O(N log N) convolution
- **Real FFT optimisation** — exploit conjugate symmetry to halve the work for real-valued inputs
- **Production use** — consider [FFTW](http://www.fftw.org/) or Intel IPP for real-time DSP workloads

---

## 📚 Further Reading

- Cooley, J.W. & Tukey, J.W. (1965). *An Algorithm for the Machine Calculation of Complex Fourier Series.* Mathematics of Computation.
- Oppenheim, A.V. & Schafer, R.W. — *Discrete-Time Signal Processing* (3rd ed.)
- [The Scientist and Engineer's Guide to DSP](http://www.dspguide.com/) — free online textbook

---

## 📄 License

Released into the public domain. Use freely for learning, research, and production.
