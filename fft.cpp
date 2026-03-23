/*
 * fft_demo.cpp
 * FFT demonstration in Digital Signal Processing using C++
 *
 * Concepts demonstrated:
 *   - Cooley-Tukey Radix-2 Decimation-In-Time (DIT) FFT
 *   - Complex arithmetic with <complex>
 *   - Frequency spectrum analysis (magnitude & phase)
 *   - Inverse FFT (IFFT)
 *   - Application: detecting dominant frequencies in a multi-tone signal
 *
 * Build: g++ -std=c++17 -O2 -o fft_demo fft_demo.cpp
 * Run:   ./fft_demo
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <string>

using Complex = std::complex<double>;
using ComplexVec = std::vector<Complex>;

const double PI = std::acos(-1.0);

// ─── Utilities ────────────────────────────────────────────────────────────────

// Check if n is a power of 2 (required for Cooley-Tukey radix-2 FFT)
bool isPowerOfTwo(size_t n) { return n > 0 && (n & (n - 1)) == 0; }

// Pad a signal with zeros up to the next power of 2
ComplexVec zeroPad(const ComplexVec& x) {
    size_t n = 1;
    while (n < x.size()) n <<= 1;
    ComplexVec padded(x);
    padded.resize(n, {0.0, 0.0});
    return padded;
}

// ─── Core FFT (Cooley-Tukey Radix-2 DIT, in-place) ───────────────────────────
//
//  The Discrete Fourier Transform (DFT) of N samples is:
//
//      X[k] = Σ_{n=0}^{N-1}  x[n] · e^{-j 2π k n / N}
//
//  Cooley-Tukey splits the N-point DFT into two N/2 DFTs (even/odd indices)
//  and combines them using "butterfly" operations, reducing complexity
//  from O(N²) to O(N log₂ N).
//
//  Butterfly:
//      E[k] + W^k · O[k]      (top output)
//      E[k] - W^k · O[k]      (bottom output)
//  where W = e^{-j 2π / N} is the twiddle factor.

void fft(ComplexVec& x, bool inverse = false) {
    size_t N = x.size();
    assert(isPowerOfTwo(N) && "FFT size must be a power of 2");

    // Bit-reversal permutation
    for (size_t i = 1, j = 0; i < N; ++i) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(x[i], x[j]);
    }

    // Cooley-Tukey butterfly stages
    for (size_t len = 2; len <= N; len <<= 1) {
        double angle = 2 * PI / len * (inverse ? 1 : -1);
        Complex wlen(std::cos(angle), std::sin(angle));   // principal root of unity

        for (size_t i = 0; i < N; i += len) {
            Complex w(1.0, 0.0);                          // current twiddle factor
            for (size_t j = 0; j < len / 2; ++j) {
                Complex u = x[i + j];
                Complex v = x[i + j + len / 2] * w;
                x[i + j]            = u + v;             // butterfly top
                x[i + j + len / 2]  = u - v;             // butterfly bottom
                w *= wlen;
            }
        }
    }

    // Normalize IFFT output by 1/N
    if (inverse)
        for (auto& val : x) val /= static_cast<double>(N);
}

// ─── Analysis Helpers ─────────────────────────────────────────────────────────

// Return magnitude spectrum (one-sided, 0..N/2)
std::vector<double> magnitudeSpectrum(const ComplexVec& X) {
    std::vector<double> mag;
    mag.reserve(X.size() / 2 + 1);
    for (size_t k = 0; k <= X.size() / 2; ++k)
        mag.push_back(std::abs(X[k]));
    return mag;
}

// Return phase spectrum in radians (one-sided)
std::vector<double> phaseSpectrum(const ComplexVec& X) {
    std::vector<double> phase;
    phase.reserve(X.size() / 2 + 1);
    for (size_t k = 0; k <= X.size() / 2; ++k)
        phase.push_back(std::arg(X[k]));
    return phase;
}

// Convert frequency bin index → actual frequency in Hz
double binToHz(size_t k, size_t N, double sampleRate) {
    return static_cast<double>(k) * sampleRate / N;
}

// Simple ASCII bar chart
void printBar(const std::string& label, double value, double maxVal, int width = 40) {
    int bars = static_cast<int>(value / maxVal * width);
    std::cout << std::setw(12) << label << " │";
    for (int i = 0; i < bars; ++i) std::cout << "█";
    std::cout << " " << std::fixed << std::setprecision(1) << value << "\n";
}

// ─── Main Demo ────────────────────────────────────────────────────────────────

int main() {
    // ── Signal parameters ──────────────────────────────────────────────────
    const double SAMPLE_RATE = 1000.0;   // Hz  (samples per second)
    const size_t N           = 256;      // must be a power of 2

    // We synthesise a signal made of three pure tones
    //   Tone A: 50 Hz,  amplitude 1.0
    //   Tone B: 120 Hz, amplitude 0.5
    //   Tone C: 200 Hz, amplitude 0.75
    const double freqA = 50.0,  ampA = 1.00;
    const double freqB = 120.0, ampB = 0.50;
    const double freqC = 200.0, ampC = 0.75;

    std::cout << "╔══════════════════════════════════════════════════╗\n";
    std::cout << "║   FFT Demo — Digital Signal Processing in C++   ║\n";
    std::cout << "╚══════════════════════════════════════════════════╝\n\n";
    std::cout << "Signal composition:\n";
    std::cout << "  Tone A: " << freqA << " Hz, amplitude " << ampA << "\n";
    std::cout << "  Tone B: " << freqB << " Hz, amplitude " << ampB << "\n";
    std::cout << "  Tone C: " << freqC << " Hz, amplitude " << ampC << "\n\n";

    // ── 1. Generate the composite time-domain signal ───────────────────────
    ComplexVec signal(N);
    for (size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / SAMPLE_RATE;
        double sample =
            ampA * std::sin(2 * PI * freqA * t) +
            ampB * std::sin(2 * PI * freqB * t) +
            ampC * std::sin(2 * PI * freqC * t);
        signal[n] = {sample, 0.0};
    }

    // Print a few samples
    std::cout << "── Time-domain samples (first 8) ──\n";
    for (size_t n = 0; n < 8; ++n) {
        std::cout << "  x[" << std::setw(2) << n << "] = "
                  << std::fixed << std::setprecision(6) << signal[n].real() << "\n";
    }

    // ── 2. Forward FFT ─────────────────────────────────────────────────────
    ComplexVec spectrum = signal;        // copy; FFT is in-place
    fft(spectrum);                       // forward transform

    auto magnitude = magnitudeSpectrum(spectrum);
    auto phase     = phaseSpectrum(spectrum);

    // Normalize magnitudes: for a real signal, single-sided peak amplitude
    // = |X[k]| * 2 / N   (factor 2 for one-sided; /N for DFT definition)
    double normFactor = 2.0 / N;

    // ── 3. Find peak frequencies ───────────────────────────────────────────
    std::cout << "\n── Frequency spectrum (peaks > 0.1 amplitude) ──\n";

    struct Peak { double freq; double amp; double phase; };
    std::vector<Peak> peaks;

    for (size_t k = 1; k < magnitude.size(); ++k) {
        double normAmp = magnitude[k] * normFactor;
        if (normAmp > 0.1) {            // threshold to ignore noise
            double freq = binToHz(k, N, SAMPLE_RATE);
            peaks.push_back({freq, normAmp, phase[k]});
        }
    }

    // Sort by amplitude descending
    std::sort(peaks.begin(), peaks.end(),
        [](const Peak& a, const Peak& b){ return a.amp > b.amp; });

    for (const auto& p : peaks) {
        std::cout << "  f = " << std::setw(7) << std::fixed << std::setprecision(2)
                  << p.freq << " Hz"
                  << "   |amplitude| = " << std::setprecision(4) << p.amp
                  << "   phase = " << std::setprecision(4) << p.phase << " rad\n";
    }

    // ── 4. ASCII magnitude spectrum bar chart ──────────────────────────────
    std::cout << "\n── Magnitude spectrum (one-sided) ──\n";

    double maxMag = *std::max_element(magnitude.begin(), magnitude.end());
    for (size_t k = 0; k < magnitude.size(); ++k) {
        double normAmp = magnitude[k] * normFactor;
        if (normAmp > 0.05) {           // only print visible lines
            double freq = binToHz(k, N, SAMPLE_RATE);
            std::string label = std::to_string(static_cast<int>(std::round(freq))) + " Hz";
            printBar(label, normAmp, 1.0, 36);
        }
    }

    // ── 5. Inverse FFT — reconstruct the original signal ──────────────────
    ComplexVec recovered = spectrum;
    fft(recovered, /*inverse=*/true);

    // Measure reconstruction error (should be essentially 0)
    double maxError = 0.0;
    for (size_t n = 0; n < N; ++n) {
        maxError = std::max(maxError, std::abs(recovered[n].real() - signal[n].real()));
    }

    std::cout << "\n── Inverse FFT (IFFT) reconstruction ──\n";
    std::cout << "  Max round-trip error: " << std::scientific
              << std::setprecision(2) << maxError
              << " (should be near machine epsilon)\n";

    std::cout << "\n── Recovered samples vs original (first 5) ──\n";
    for (size_t n = 0; n < 5; ++n) {
        std::cout << "  x[" << n << "] original = " << std::fixed << std::setprecision(8)
                  << signal[n].real()
                  << "   recovered = " << recovered[n].real() << "\n";
    }

    // ── 6. Application: simple low-pass filter in frequency domain ─────────
    std::cout << "\n── Low-pass filter (cutoff = 100 Hz) ──\n";
    ComplexVec filtered = spectrum;
    const double cutoff = 100.0;

    // Zero out bins above the cutoff (rectangular window LPF)
    for (size_t k = 0; k < filtered.size(); ++k) {
        double freq = binToHz(k < filtered.size() / 2 ? k : filtered.size() - k,
                              N, SAMPLE_RATE);
        if (freq > cutoff) filtered[k] = {0.0, 0.0};
    }

    // IFFT back to time domain
    fft(filtered, /*inverse=*/true);

    std::cout << "  Frequencies above " << cutoff << " Hz removed.\n";
    std::cout << "  Filtered sample x[10] = " << std::fixed << std::setprecision(6)
              << filtered[10].real() << "\n";

    // The only surviving tone should be freqA (50 Hz)
    // Re-run FFT on filtered signal and show peak
    ComplexVec filtSpec = filtered;
    fft(filtSpec);
    auto filtMag = magnitudeSpectrum(filtSpec);
    double filtMax = *std::max_element(filtMag.begin(), filtMag.end());
    size_t peakBin = std::max_element(filtMag.begin(), filtMag.end()) - filtMag.begin();
    std::cout << "  Dominant frequency after filtering: "
              << std::fixed << std::setprecision(1)
              << binToHz(peakBin, N, SAMPLE_RATE) << " Hz"
              << " (expected " << freqA << " Hz)\n";

    std::cout << "\n✓ Demo complete.\n";
    return 0;
}