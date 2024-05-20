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

I had a few published research regards to high speed implementation of NTT. 
For short, 

1. In table 4 of my paper[2], a vectorize version of NTT faster than the C version with  compiler `-O3` (highest optimization) by `6 to 8` times. 
2. NTT has the best time complexity `O(nlogn)`. 
3. Compiler does not understand how to vectorize NTT properly. 

How can 1 it impact the performance of ZKVM? 

1. Batch processing SIMD NTT instructions -> Lead to at least 4x faster. 
2. Higher level of parallelism, the current implementation does not make all my CPU cores full. There is not enough task for the CPU to run! If this is done properly, expect to be at least 2x faster. 

With these estimation, I expect there is room for 8x improvement in speed. The prove time from 107 minutes can be reduced down to 13 minutes. It's 94 minutes different!. 

## How to run

`cargo nexus verify`.



[1] Huang, Junhao, et al. "Improved Plantard arithmetic for lattice-based cryptography." IACR Transactions on Cryptographic Hardware and Embedded Systems 2022.4 (2022): 614-636.
[2] "Fast Falcon Signature Generation and Verification Using ARMv8 NEON Instructions" https://csrc.nist.gov/csrc/media/Events/2022/fourth-pqc-standardization-conference/documents/papers/fast-falcon-signature-generation-and-verification-pqc2022.pdf 