# PrimeLister (v3.1.1 – v3.4.1)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24451603.svg)](https://doi.org/10.5281/zenodo.24451603)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![SIMD](https://img.shields.io/badge/SIMD-AVX2%20%7C%20AVX--512-orange.svg)
![OpenMP](https://img.shields.io/badge/Parallel-OpenMP-green.svg)

**PrimeLister** is a single-header/single-file C++17 parallel segmented prime sieve designed for **full 64-bit RAM materialization** (generating and storing every prime as a `uint64_t` in memory, not just counting).

All versions from `v3.1.1` through `v3.4.1` strictly preserve the author's original **modulo-6 candidate mapping and composite-marking equations**, demonstrating how far cache residency, branchless SIMD extraction (`VPCOMPRESSQ`), and non-temporal DRAM streaming stores (`_mm512_stream_si512` / `_mm256_stream_si256`) can push a compact Mod-6 sieve on modern CPUs.

> **Full Technical Report (Rev. 2.1) & Raw Audit Data:** [Zenodo Record 24451603](https://zenodo.org/records/24451603)

---

## Key Highlights (AMD Ryzen 9 9955HX, 16C/32T, DDR5-5600, MSVC AVX-512)

* **32-Thread RAM Materialization vs. `primesieve 12.18`:** PrimeLister `v3.3`, `v3.4`, and `v3.4.1` outperform an OpenMP-parallelized `primesieve::generate_primes()` harness across all tested limits ($N = 10^9$ to $1.5 \times 10^{10}$) by **+22% to +38% mean** (**+19% to +51% median**), and outperform the out-of-the-box single-threaded `primesieve::generate_primes()` API by **7.3× to 9.1×**.
* **Up to 4.13 billion primes/s materialized to RAM:** `v3.4` reaches **3.712 B/s mean** (3.999 B/s median, 4.132 B/s max) at $N = 10^9$ and **3.557 B/s** at $N = 10^{10}$.
* **5.36 GB of primes in 0.201 s ($N = 1.5 \times 10^{10}$):** `v3.4.1` materializes all **670,180,516 primes** up to 15 billion in **0.2015 s** (**3.332 B/s**, 26.6 GB/s net payload write rate), beating parallelized `primesieve 12.18` by **+29.0% mean (+36.4% median)**.
* **Honest trade-off:** In **32-thread counting-only mode** (`--count`, no RAM output) and on **large single-threaded ranges**, Kim Walisch's `primesieve 12.18` remains faster (~7.0–7.2 B/s vs. ~5.1–5.8 B/s counting) due to its modulo-210 wheel and AVX-512 `VPOPCNTDQ`.

---

## Preserved Mathematical Core

After handling `2` and `3` explicitly, PrimeLister maps all integers $\ge 5$ coprime to $6$ ($6k \pm 1$) to zero-based indices:

```cpp
candidate(i) = 3 * i + ((i & 1) ? 4 : 5);
```

For each base prime $p \le \lfloor\sqrt{N}\rfloor$, composites are marked via exactly two progressions of stride `2 * p`:

1. **Square progression:** `start = (p * p) / 3 - 1`, `step = 2 * p`
2. **Cofactor progression:** `start = (p * other) / 3 - 1`, `step = 2 * p`, where `other = p + (p % 6 == 5 ? 2 : 4)`

---

## Version Evolution

| Version | Segment | Presieving | Inner Loops & Bucket Sieve | Extraction & Stores |
| :--- | :---: | :--- | :--- | :--- |
| **`v3.1.1`** | 128 KiB | None | Generic 64-bit word loop | Scalar `tzcnt`, `std::vector::resize`, SWAR popcnt |
| **`v3.2`** | 32 KiB | 8 patterns (`5..73`) | Fused 8-mark/branch + bucket ($\ge 131{,}072$) | `CANDIDATE_LUT`, uninitialized chunks, HW `POPCNT` |
| **`v3.3`** | 32 KiB | 16 tiled AVX2 (`5..163`) | Switch-free 8-lane stride-$p$ + flat word + bucket | 64 B aligned chunks + 4 KiB L1 stage + **AVX2 NT stores** |
| **`v3.4`** | 24 KiB | 16 tiled AVX-512 `VPTERNLOGQ` | Switch-free 8-lane + 2× flat + bucket ($\ge 98{,}304$) | Branchless **AVX-512 `VPCOMPRESSQ`** + **64 B NT stores** |
| **`v3.4.1`** | 32 KiB | 16 tiled AVX-512 `VPTERNLOGQ` | Switch-free 8-lane + flat + shift-bucket ($\ge 131{,}072$) | Branchless **AVX-512 `VPCOMPRESSQ`** + **64 B NT stores** |

*(Full AVX2 and portable scalar fallbacks are included automatically when compiled without `/arch:AVX512` or `-mavx512f`.)*

---

## Matched Head-to-Head Benchmark vs. `primesieve 12.18`

Measured in a single executable (`HeadToHead Benchmarker v2.1`, isolated sessions per $N$, $n = 5$ interleaved randomized rounds, pre-flight prime-by-prime check at $N = 10^6$, Windows 11, MSVC `195136257`, `/O2 /Oi /Ot /GL /arch:AVX512 /openmp:llvm`).

### 1. 32-Thread Full 64-Bit RAM Materialization (Mean Rate, Billion Primes/s)

| Limit $N$ (Payload) | `PL 3.4.1` (32K) | `PL 3.4` (24K) | `PL 3.3` (AVX2) | `primesieve 12.18` (32T OMP) | `primesieve` (1T API) | Best PL vs. `primesieve` (32T) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`1e9`** (0.41 GB) | 3.282 | **3.712** *(Med 3.999)* | 3.656 | 2.688 | 0.408 | **+38.1%** *(Med +50.6%)* |
| **`5e9`** (1.88 GB) | 3.324 | **3.523** | 3.250 | 2.660 | 0.440 | **+32.4%** |
| **`10e9`** (3.64 GB) | 3.474 | **3.557** *(Med 3.625)* | 3.445 | 2.622 | 0.461 | **+35.7%** *(Med +41.7%)* |
| **`15e9`** (5.36 GB) | **3.332** *(0.2015 s)* | 2.957 | 3.267 | 2.582 | 0.456 | **+29.0%** *(Med +36.4%)* |

### 2. Single-Thread RAM (1T) & Pure Counting (`--count`, 32T)

| Limit $N$ | `PL 3.4.1` (1T RAM) | `primesieve` (1T RAM) | `PL 3.4.1` (32T Count) | `PL 3.4` (32T Count) | `primesieve` (32T Count) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`1e9`** | 0.396 B/s | **0.408 B/s** | 5.151 B/s | 4.810 B/s | **5.967 B/s** |
| **`5e9`** | 0.295 B/s | **0.440 B/s** | 5.663 B/s | 5.785 B/s | **7.163 B/s** |
| **`10e9`** | 0.304 B/s | **0.461 B/s** | 5.311 B/s | 5.479 B/s | **7.006 B/s** |
| **`15e9`** | 0.285 B/s | **0.456 B/s** | 5.174 B/s | 4.335 B/s | **7.083 B/s** |

**Why PrimeLister wins 32T RAM materialization:** Standard `std::vector` stores trigger DRAM *Write-Allocate* (Read-For-Ownership) traffic, capping `primesieve_OMP_MT` at ~20.5–21.4 GB/s net payload write throughput. PrimeLister's 4 KiB L1 staging buffer combined with 64-byte non-temporal streaming stores bypasses cache write-allocation and sustains **26.6–29.2 GB/s** net payload into DDR5-5600.

**Why `v3.4` (24 KiB) and `v3.4.1` (32 KiB) cross over at $1.5 \times 10^{10}$:** Up to $10^{10}$, `v3.4`'s 24 KiB segment fits both SMT threads cleanly into Zen 5's 48 KiB L1d cache ($2 \times 24 = 48$ KiB). At $N = 1.5 \times 10^{10}$ ($\sqrt{N} = 122{,}474$), `v3.4` crosses its bucket threshold (`98,304`) and incurs non-power-of-two segment division, whereas `v3.4.1`'s power-of-two 32 KiB segment (`BUCKET_LIMIT = 131,072`) stays division-free and gains +12.7% in 32T RAM (+19.3% in counting).

---

## Build & Usage

### MSVC (Visual Studio 2022/2026, x64 Native Tools Command Prompt)
```bat
cl /O2 /Oi /Ot /GL /arch:AVX512 /openmp:llvm /std:c++17 /EHsc PrimeLister_v3.4.1.cpp /link /LTCG
```
*(Use `/arch:AVX2` on CPUs without AVX-512.)*

### GCC / Clang (Linux / MSYS2)
```bash
g++ -O3 -march=native -fopenmp -std=c++17 PrimeLister_v3.4.1.cpp -o primelister
```

### CLI Options
```bash
# Compute and materialize all primes up to 10 billion (auto-calibrated threads), no disk output
./primelister 10000000000 --no-output

# Fixed 32 threads, write tab-separated table (1000 columns) to PrimeList.txt
./primelister 1000000000 --threads 32 --cols 1000 --output PrimeList.txt

# Counting-only mode (no RAM list allocated)
./primelister 15000000000 --threads 32 --count

# Disable non-temporal streaming stores (for A/B bandwidth experiments)
./primelister 10000000000 --threads 32 --no-output --no-nt
```

---

## Publications & Provenance

* **PrimeLister 3.2–3.4.1 Technical Report (Rev. 2.1, Nov 2026):** [https://zenodo.org/records/24451603](https://zenodo.org/records/24451603)
* **PrimeLister 3.1.1 Technical Report (Rev. 1.1, Sep 2026):** [https://zenodo.org/records/22979590](https://zenodo.org/records/22979590)
* **PrimeLister 2.9 Publication (Jan 2026):** [https://zenodo.org/records/18428535](https://zenodo.org/records/18428535)

## Authorship & License

Created by **Tobias Jung** (`info@scyrus.de`). Generative AI assisted with C++ SIMD/cache optimization, benchmark harness engineering, raw-log auditing, and documentation under author direction and verification.

Released under **CC0 1.0 Universal** (or adjust to match your Zenodo deposit license).