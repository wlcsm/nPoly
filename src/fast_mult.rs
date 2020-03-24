use crate::algebras::*;
use crate::fft::*;
use crate::mathutils::*;
use crate::polyu::*;
use crate::algebras::complex::*;
use crate::algebras::integers::*;

pub trait FastMult: Ring {
    fn fast_mult(&self, b: &Self, m: usize) -> Self;
}

// TODO these two should also be a macro
fn to_coeffs_complex(input: &Monomials<ZZ>, n: usize) -> Vec<CC> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<CC> = Vec::with_capacity(n);
    for (coeff, deg) in input.iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(deg, CC::zero());
        result.push(CC::from_re(coeff.0));
    }
    // Pad the rest
    result.resize(n, CC::zero());
    result
}

fn to_coeffs<T: SupportsFFT>(input: &Monomials<T>, n: usize) -> Vec<T> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<T> = Vec::with_capacity(n);
    for (coeff, deg) in input.iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(deg, <T>::zero());
        result.push(coeff);
    }
    // Pad the rest
    result.resize(n, <T>::zero());
    result
}

// TODO these two implementations should be one macro
impl FastMult for PolyU<CC> {

    fn fast_mult(&self, other: &Self, m: usize) -> Self {

        let n = next_npow(self.deg() + other.deg() + 1, m);
        let mut a_sig = to_coeffs(&self.terms, n);
        let mut b_sig = to_coeffs(&other.terms, n);

        // Infix on a_sig
        eval_interp(&mut a_sig[..], &mut b_sig[..], m).unwrap();

        // Need to normalise it here
        for x in a_sig.iter_mut() {
            x.0 /= n as f64
        }

        // Convert back into polynomial type
        PolyU::from_coeff(None, a_sig).unwrap()
    }
}

impl FastMult for PolyU<ZZ> {

    fn fast_mult(&self, other: &Self, m: usize) -> Self {

        let n = next_npow(self.deg() + other.deg() + 1, m);
        let mut a_sig = to_coeffs_complex(&self.terms, n);
        let mut b_sig = to_coeffs_complex(&other.terms, n);

        eval_interp(&mut a_sig[..], &mut b_sig[..], m).unwrap();

        // Because we converted it to complex for the ROU
        // We also need to normalise it here
        let c_parsed = a_sig.into_iter()
                            .map(|x| 
                                ZZ((x.0 / n as f64).re.round() as i32)
                            ).collect();

        // Convert back into polynomial type
        PolyU::from_coeff(None, c_parsed).unwrap()
    }
}