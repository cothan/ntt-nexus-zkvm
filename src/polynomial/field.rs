use core::ops::{Add, Mul, Sub};

pub type Integer = i16;

/// Represents an element of GF(q). Even though `q` is 16 bits wide, we use a wider uint type so that we can defer modular reductions.
#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd, Eq, Ord)]
pub struct FieldElement(pub Integer);

impl FieldElement {
    // Fermat prime, special modular reduction
    pub const Q: Integer = 257;
    pub const Q32: i32 = Self::Q as i32;
    pub const Q64: i64 = Self::Q as i64;
    pub const QINV: i64 = 4278255361;

    /**
     * Reduces the input value using Fermat special modular reduction.
     *
     * # Arguments
     * * `x` - The integer to be reduced.
     *
     * # Returns
     * The reduced integer.
     */
    pub fn small_reduce(x: Integer) -> Integer {
        let a = x >> 8;
        let b = x & 0xff;
        let r = b - a;
        if r < 0 {
            r + Self::Q
        } else {
            r
        }
    }

    /**
     * Precomputes a constant for Plantard multiplication method.
     *
     * # Arguments
     * * `b` - The integer for which the constant is precomputed.
     *
     * # Returns
     * The precomputed constant as a 32-bit integer.
     */
    fn _plantard_precompute_constant(b: Integer) -> i32 {
        let mut b_prime = i64::from(b);
        b_prime <<= 32;
        b_prime = -b_prime;
        b_prime %= Self::Q64;

        b_prime *= Self::QINV;
        b_prime &= 0xffffffff;

        b_prime as i32
    }

    /**
     * Multiplies an integer by a precomputed constant using Plantard method.
     *
     * # Arguments
     * * `a` - The integer to be multiplied.
     * * `bq_prime` - The precomputed constant.
     *
     * # Returns
     * The result of the multiplication as an integer.
     */
    fn plantard_mul_with_constant(a: Integer, bq_prime: i32) -> Integer {
        let mut abq_prime: i32 = bq_prime.wrapping_mul(a as i32);

        abq_prime >>= 16;
        abq_prime += 1 << 6;
        abq_prime = abq_prime * Self::Q32;
        abq_prime >>= 16;

        abq_prime.try_into().unwrap()
    }

    /**
     * Multiplies two integers using Barrett reduction.
     * Arguments
     * `a` - The first integer.
     * `b` - The second integer.
     *
     * Returns
     * The result of the multiplication as an integer.
     */
    fn barrett_mul(a: Integer, b: Integer) -> Integer {
        let c = (a as i32).wrapping_mul(b as i32) % Self::Q32;

        c.try_into().unwrap()
    }

    /**
     * Creates a `FieldElement` with a value of zero.
     *
     * # Returns
     * A `FieldElement` set to zero.
     */
    pub fn zero() -> Self {
        Self(0)
    }
}

impl Add<FieldElement> for FieldElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(Self::small_reduce(self.0 + rhs.0))
    }
}

impl Sub<FieldElement> for FieldElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(Self::small_reduce(self.0 - rhs.0))
    }
}

impl Mul<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn mul(self, rhs: FieldElement) -> Self {
        Self(Self::barrett_mul(self.0, rhs.0))
    }
}

impl Mul<i32> for FieldElement {
    type Output = FieldElement;

    fn mul(self, rhs: i32) -> Self {
        Self(Self::plantard_mul_with_constant(self.0, rhs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = FieldElement(123);
        let b = FieldElement(197);
        let result = a + b;
        assert_eq!(result, FieldElement(63));
    }

    #[test]
    fn test_sub() {
        let a = FieldElement(123);
        let b = FieldElement(197);
        let result = a - b;
        assert_eq!(result, FieldElement(183));
    }

    #[test]
    fn test_mul() {
        let a = FieldElement(123);
        let b = FieldElement(197);
        let result = a * b;
        assert_eq!(result, FieldElement(73));
    }

    #[test]
    fn test_reduce() {
        let x = 513;
        let result = FieldElement::small_reduce(x);
        assert_eq!(result, 256);
    }
}
