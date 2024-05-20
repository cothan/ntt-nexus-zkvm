pub mod field;

use core::ops::{Add, Mul, Sub};
use field::{FieldElement, Integer};

#[derive(PartialEq, Debug, Eq, PartialOrd, Ord)]
/**
 * Represents a polynomial with coefficients in a finite field.
 *
 * Fields:
 * - `coeffs`: An array of `FieldElement`s representing the coefficients of the polynomial.
 *             The length of the array is fixed as `Polynomial::LENGTH`.
 */
pub struct Polynomial {
    pub coeffs: [FieldElement; Polynomial::LENGTH],
}
impl FromIterator<FieldElement> for [FieldElement; Polynomial::LENGTH] {
    /**
     * Implements the `FromIterator` trait for converting an iterator of `FieldElement`s
     * into an array of `FieldElement`s of fixed size `Polynomial::LENGTH`.
     *
     * Parameters:
     * - `iter`: An iterator that yields `FieldElement` items.
     *
     * Returns:
     * `[FieldElement; Polynomial::LENGTH]` - An array containing the `FieldElement` items
     * from the iterator, with any remaining slots filled with the zero `FieldElement`.
     */
    fn from_iter<I: IntoIterator<Item = FieldElement>>(iter: I) -> Self {
        let mut result = [FieldElement::zero(); Polynomial::LENGTH]; // Assuming zero is a valid default value
        for (i, element) in iter.into_iter().enumerate() {
            result[i] = element;
        }
        result
    }
}

impl Add<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    /**
     * Adds two polynomials and returns the resulting polynomial.
     *
     * Parameters:
     * - `self`: Reference to the first polynomial.
     * - `rhs`: Reference to the second polynomial to be added to the first.
     *
     * Returns:
     * `Polynomial` - The polynomial resulting from the addition of `self` and `rhs`.
     */
    fn add(self, rhs: &Polynomial) -> Polynomial {
        Polynomial {
            coeffs: self
                .coeffs
                .iter()
                .zip(rhs.coeffs.iter())
                .map(|(&x, &y)| x + y)
                .collect(),
        }
    }
}

impl Sub<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    /**
     * Subtracts the second polynomial from the first and returns the resulting polynomial.
     *
     * Parameters:
     * - `self`: Reference to the first polynomial.
     * - `rhs`: Reference to the second polynomial to be subtracted from the first.
     *
     * Returns:
     * `Polynomial` - The polynomial resulting from the subtraction of `rhs` from `self`.
     */
    fn sub(self, rhs: &Polynomial) -> Polynomial {
        Polynomial {
            coeffs: self
                .coeffs
                .iter()
                .zip(rhs.coeffs.iter())
                .map(|(&x, &y)| x - y)
                .collect(),
        }
    }
}

impl Mul<&Polynomial> for FieldElement {
    type Output = Polynomial;

    /**
     * Multiplies a polynomial by a field element and returns the resulting polynomial.
     *
     * Parameters:
     * - `self`: The field element by which the polynomial is to be multiplied.
     * - `rhs`: Reference to the polynomial to be multiplied.
     *
     * Returns:
     * `Polynomial` - The polynomial resulting from the multiplication of `self` and `rhs`.
     */
    fn mul(self, rhs: &Polynomial) -> Polynomial {
        Polynomial {
            coeffs: rhs.coeffs.iter().map(|&x| self * x).collect(),
        }
    }
}

/**
 * Precompute the twiddle factors, this help accelerate the performance of Number Theoretic Transform
 */
const ZETA_POW_BITREV: [i32; Polynomial::LENGTH] = {
    const ZETA: Integer = Polynomial::ROOT_OF_UNITY;
    const fn bitrev6(x: usize) -> usize {
        ((x >> 5) & 1)
            | (((x >> 4) & 1) << 1)
            | (((x >> 3) & 1) << 2)
            | (((x >> 2) & 1) << 3)
            | (((x >> 1) & 1) << 4)
            | ((x & 1) << 5)
    }

    // Compute the powers of zeta
    let mut pow = [0i64; Polynomial::LENGTH];
    let mut i = 0;
    let mut curr = 1i32;
    while i < Polynomial::LENGTH {
        pow[i] = curr as i64;
        i += 1;
        curr = curr * (ZETA as i32) % FieldElement::Q32;
    }

    let mut pow_bitrev = [0i32; Polynomial::LENGTH];
    let mut i = 0;
    while i < Polynomial::LENGTH {
        let mut b_prime = pow[bitrev6(i)];
        b_prime <<= 32;
        b_prime = -b_prime;
        b_prime %= FieldElement::Q64;

        b_prime *= FieldElement::QINV;
        b_prime &= 0xffffffff;

        pow_bitrev[i] = b_prime as i32;
        i += 1;
    }
    pow_bitrev
};
impl Polynomial {
    pub const LENGTH: usize = 64;
    pub const ROOT_OF_UNITY: Integer = 9;
    /**
     * Performs the Number Theoretic Transform (NTT) on the polynomial.
     *
     * This method modifies the polynomial's coefficients in-place using the NTT algorithm.
     * The NTT is similar to the Fast Fourier Transform (FFT) but operates in finite fields.
     *
     * Returns:
     * `NttPolynomial` - The transformed polynomial in the NTT domain.
     */
    pub fn ntt(&self) -> NttPolynomial {
        let mut k: usize = 1;

        let mut f = self.coeffs;

        for len in [32, 16, 8, 4, 2, 1] {
            for start in (0..Self::LENGTH).step_by(2 * len) {
                let zeta: i32 = ZETA_POW_BITREV[k];
                k += 1;

                for j in start..(start + len) {
                    let t = f[j + len] * zeta;
                    f[j + len] = f[j] - t;
                    f[j] = f[j] + t;
                }
            }
        }

        f.into()
    }
}

#[derive(PartialEq, Debug, Eq, PartialOrd, Ord)]
/**
 * Represents a polynomial in the Number Theoretic Transform (NTT) domain.
 *
 * Fields:
 * - `coeffs`: An array of `FieldElement`s representing the coefficients of the polynomial.
 *             The length of the array is fixed as `Polynomial::LENGTH`.
 */
pub struct NttPolynomial {
    pub coeffs: [FieldElement; Polynomial::LENGTH],
}

impl From<[FieldElement; Polynomial::LENGTH]> for NttPolynomial {
    fn from(arr: [FieldElement; Polynomial::LENGTH]) -> Self {
        NttPolynomial { coeffs: arr }
    }
}

impl NttPolynomial {
    pub fn init() -> Self {
        Self {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        }
    }
    /**
     * This function performs the inverse Number Theoretic Transform (NTT) on the polynomial.
     * It modifies the NTT polynomial's coefficients using the inverse NTT algorithm.
     *
     * Returns:
     * `Polynomial` - The transformed polynomial back in the normal coefficient representation.
     */
    pub fn ntt_inverse(&mut self) -> Polynomial {
        let mut f = self.coeffs;

        let mut k = ZETA_POW_BITREV.len() - 1;
        for len in [1, 2, 4, 8, 16, 32] {
            for start in (0..Polynomial::LENGTH).step_by(2 * len) {
                let zeta = ZETA_POW_BITREV[k];
                k -= 1;

                for j in start..(start + len) {
                    let t = f[j];
                    f[j] = t + f[j + len];
                    f[j + len] = (f[j + len] - t) * zeta;
                }
            }
        }

        for i in 0..f.len() {
            f[i] = f[i] * FieldElement(253);
        }

        Polynomial { coeffs: f }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_zero() {
        let a = Polynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };
        let b = Polynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };

        let c = a.ntt().ntt_inverse();

        assert_eq!(b, c);
    }

    #[test]
    fn test_ntt_then_inverse_ntt() {
        let mut a = Polynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };
        let mut b = Polynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };

        for i in 0..Polynomial::LENGTH {
            a.coeffs[i] = FieldElement(i as i16);
            b.coeffs[i] = FieldElement(i as i16);
        }

        let mut c = a.ntt().ntt_inverse();

        for i in 0..Polynomial::LENGTH {
            c.coeffs[i] = FieldElement(FieldElement::small_reduce(c.coeffs[i].0));
        }

        assert_eq!(b, c);
    }

    #[test]
    fn test_inverse_ntt_then_ntt() {
        let mut a = NttPolynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };
        let mut b = NttPolynomial {
            coeffs: [FieldElement::zero(); Polynomial::LENGTH],
        };

        for i in 0..Polynomial::LENGTH {
            a.coeffs[i] = FieldElement(i as i16);
            b.coeffs[i] = FieldElement(i as i16);
        }

        let c = a.ntt_inverse().ntt();

        assert_eq!(b, c);
    }
}
