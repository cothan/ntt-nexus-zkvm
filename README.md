# Number Theoretic Transform in Nexus ZKVM

## Original idea

I plan to implement Post-quantum signature algorithm to verify a signature. 
However, given that ZKVM takes so much time to computer NTT size `n = 16`, the wait is too long so I need to scale down the work. 
Real world PQC digital signature is `n = 512, 1024`, SHA3 hashing, and lattice-based polynomial arithmetic.

## Number Theoretic Transform (NTT)

The Number Theoretic Transform (NTT) is a mathematical technique used for efficiently performing polynomial multiplication in finite fields. It is the discrete Fourier transform (DFT) adapted to arithmetic in a finite field. The NTT is particularly useful in cryptographic applications, such as in lattice-based cryptography and homomorphic encryption.

The NTT is a powerful algorithm for efficient polynomial multiplication in finite fields, leveraging properties of roots of unity and modular arithmetic to achieve significant performance improvements in computational tasks.

## Plantard Multiplication

This work use the new state of the art multiplication with constant: Plantard Multiplication[1].
The detail can be found in the improved Plantard[1]. 

I implement it here to speed up the implementation of NTT on Nexus ZKVM: reduce computation, hence faster. 

## Performance 

- Prove: 1 hour 47 minutes = 107 minutes.
- Verify: 5 seconds.

![image](./Screenshot%202024-05-19%20at%2011.16.16 PM.png)

## Real world performance of NTT

I have published several research papers on high-speed implementations of the Number Theoretic Transform (NTT).

In brief:

- In Table 4 of my paper [2], a vectorized version of NTT outperforms the C version compiled with the highest optimization level (-O3) by a factor of 6 to 8.
- NTT achieves an optimal time complexity of O(nlog⁡n)O(nlogn).
- Current compilers do not adequately vectorize NTT. A manual implementation of NTT significantly surpasses the performance of a pure C implementation. Consequently, one should not evaluate the performance of NTT solely based on compiler optimization.

Potential Impact on ZKVM Performance:

- Implementing batch processing with SIMD NTT instructions could result in at least a 4x speedup.
- Increasing the level of parallelism is crucial. The current implementation does not fully utilize all CPU cores, as there are insufficient tasks for the CPU to execute. Proper parallelism could yield at least a 2x improvement in performance.

With these estimations, there is potential for an 8x improvement in speed. This enhancement could reduce the proof time from 107 minutes to approximately 13 minutes, a reduction of 94 minutes.

## How to run

`cargo nexus verify`.

## References

[1] Huang, Junhao, et al. "Improved Plantard arithmetic for lattice-based cryptography." IACR Transactions on Cryptographic Hardware and Embedded Systems 2022.4 (2022): 614-636.

[2] "Fast Falcon Signature Generation and Verification Using ARMv8 NEON Instructions" https://csrc.nist.gov/csrc/media/Events/2022/fourth-pqc-standardization-conference/documents/papers/fast-falcon-signature-generation-and-verification-pqc2022.pdf 