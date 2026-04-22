# BCH(63, 36, 5) Error Correction

A lightweight, self-contained C implementation of a **BCH(63, 36, 5)** error-correcting code over GF(2⁶). It can encode 36-bit messages into 63-bit codewords and correct up to **5 random bit errors** during decoding.

---

## What is BCH?

BCH (Bose–Chaudhuri–Hocquenghem) codes are a family of cyclic error-correcting codes. This implementation uses the following parameters:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `n` | 63 | Total codeword length (bits) |
| `k` | 36 | Message length (bits) |
| `t` | 5 | Max correctable bit errors |
| Parity bits | 27 | Redundancy added for error correction |
| Field | GF(2⁶) | Galois Field with 63 non-zero elements |
| Primitive poly | `x⁶ + x + 1` (`0x43`) | Defines the field arithmetic |
| Generator poly | `0x86e8113` (degree 27) | Defines the code |

---

## How It Works

The encode–transmit–decode pipeline:

```
 36-bit data
     │
     ▼
 ┌─────────┐      63-bit codeword
 │  Encode │ ──────────────────────────► [transmission / storage]
 └─────────┘                                       │
                                                   │ (up to 5 bit flips)
                                                   ▼
                                        ┌──────────────────────┐
                                        │  Compute Syndromes   │
                                        └──────────┬───────────┘
                                                   │
                                        ┌──────────▼───────────┐
                                        │  Error Locator (BMA) │
                                        └──────────┬───────────┘
                                                   │
                                        ┌──────────▼───────────┐
                                        │  Correct Codeword    │
                                        └──────────┬───────────┘
                                                   │
                                                   ▼
                                              36-bit data  ✓
```

**Decoding uses the Berlekamp–Massey Algorithm (BMA)** to find the error-locator polynomial σ(x), then evaluates it at every field element to locate and flip the erroneous bits.

---

## Project Structure

```
.
├── BCH.h        # Constants, GF parameters, and function declarations
├── BCH.c        # Full BCH implementation (GF arithmetic, encode, decode)
└── testing.c    # Test driver — injects 5 errors and verifies correction
```

---

## Build & Run

**Requirements:** Any C99-compatible compiler (GCC, Clang, MSVC).

```bash
# Compile
gcc -o test testing.c BCH.c -std=c99

# Run
./test
```

**Expected output:**

```
original data  = 0x34fd8a0c
codeword       = 0x4f3d5c34fd8a0c...
error codeword = ...            ← 5 bits flipped
syndromes      = xx xx xx xx xx xx xx xx xx xx
sigma          = 01 xx xx xx xx xx
fix            = ...            ← corrected
expected       = ...
match          = YES
og_data        = 0x34fd8a0c
expected       = 0x34fd8a0c
```

---

## API Reference

```c
// 1. Must be called once before anything else
void bch_init(void);

// 2. Encode a 36-bit message into a 63-bit codeword
uint64_t bch_encode(uint32_t data);

// 3. Compute the 2t syndromes of a (possibly corrupted) codeword
void bch_compute_syndromes(uint64_t codeword, uint8_t syndromes[2 * BCH_T]);

// 4. Run BMA to find the error-locator polynomial sigma
void bch_syn_error_locator(const uint8_t syndromes[2 * BCH_T], uint8_t sigma[BCH_T + 1]);

// 5. Evaluate sigma at each field element and flip error bits
uint64_t bch_syn_correction(uint64_t codeword, const uint8_t sigma[BCH_T + 1]);

// 6. Extract the original 36-bit message from a corrected codeword
uint32_t bch_original_data(uint64_t corrected_codeword);
```

---

## Limitations

- Corrects **up to 5** bit errors. More than 5 errors will produce a wrong result silently — no failure detection is built in.
- Designed for **GF(2⁶) / BCH(63,36,5)** specifically. Changing `n`, `k`, or `t` requires a new generator polynomial and field tables.
- Not optimized for throughput; intended as a clear reference implementation.

---


