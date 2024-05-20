#![no_std]
#![no_main]

mod polynomial;

use polynomial::{NttPolynomial, Polynomial};

#[nexus_rt::main]
fn main() {
    let mut a = NttPolynomial::init();
    let mut b = NttPolynomial::init();

    for i in 0..Polynomial::LENGTH {
        a.coeffs[i] = crate::polynomial::field::FieldElement(i as i16);
        b.coeffs[i] = crate::polynomial::field::FieldElement(i as i16);
    }

    let c = a.ntt_inverse().ntt();

    assert_eq!(b, c);
}
