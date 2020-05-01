use crate::algebras::*;
use crate::fft::*;
use crate::polyu::*;
use crate::algebras::polyring::*;
use crate::mathutils::log2_unchecked;

pub trait FastMult {
    fn fast_mult(&self, b: &Self) -> Self;
}

// // TODO these two should also be a macro
// fn to_coeffs_complex(input: &Vec<Term<ZZ>>, n: usize) -> Vec<CC> {
//     // Expands the input into the expanded coefficient vector (coerced into complex)
//     // Then padded with zeros to length n

//     let mut result: Vec<CC> = Vec::with_capacity(n);
//     for Term { coeff, deg } in input {
//         // Fill the gap between monomials with zeros, then add the monomial
//         result.resize(deg.0 as usize, CC::zero());
//         result.push(CC::from_re(coeff.0));
//     }
//     // Pad the rest
//     result.resize(n, CC::zero());
//     result
// }

fn to_coeffs<P: PolyRing>(input: &Vec<Term<P>>, n: usize) -> Vec<P::Coeff> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<P::Coeff> = Vec::with_capacity(n);
    for Term { coeff, deg } in input.iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(deg.tot_deg(), <P::Coeff>::zero());
        result.push(*coeff);
    }
    // Pad the rest
    result.resize(n, <P::Coeff>::zero());
    result
}

// TODO these two implementations should be one macro
impl<'a, T: SupportsFFT> FastMult for PolyU<'a, T> {

    fn fast_mult(&self, other: &Self) -> Self {

        let n = (self.deg() + other.deg() + 1).next_power_of_two();
        let mut a_sig = to_coeffs(&self.terms, n);
        let mut b_sig = to_coeffs(&other.terms, n);

        // Infix on a_sig
        eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

        // Need to normalise it here
        for x in a_sig.iter_mut() {
            x.divby2(log2_unchecked(n))
        }

        // Convert back into polynomial type
        Poly::from_coeff(self.ring, a_sig)
    }
}

// Frozen this at the moment becaue it requires a change of type from CC to ZZ
// which I don't feel like doing at the moment

// impl<'a, P: PolyRing<Coeff=CC, Var=Univariate>> FastMult for Poly<'a, P> {

//     fn fast_mult(&self, other: &Self) -> Self {

//         let n = (self.deg() + other.deg() + 1).next_power_of_two();
//         let mut a_sig = to_coeffs_complex(&self.terms, n);
//         let mut b_sig = to_coeffs_complex(&other.terms, n);

//         eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

//         // Because we converted it to complex for the ROU
//         // We also need to normalise it here
//         let c_parsed = a_sig.into_iter()
//                             .map(|x| 
//                                 ZZ((x.0 / n as f64).re.round() as i32)
//                             ).collect();

//         // Convert back into polynomial type
//         Poly::from_coeff(self.ring, c_parsed)
//     }
// }