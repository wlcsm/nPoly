use crate::polym::*;
use crate::polyu::*;
use num_traits::Zero;

/// This will be used to hold the data to perform Kronecker substitution
/// The i^th element corresponds to the value that was substituted into the i^th variable
#[derive(Clone)]
pub struct KroneckerData<'a, P: PolyRing> {
    sub_values: GenericArray<usize, P::NumVar>,
    ring: &'a P,
}

pub const FIRST_100_PRIMES: [usize; 100] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
    401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541,
];

/// Uses Kronecker substitution to convert the multivariate polynomial to a univariate one
///
/// Need to find the maximum of the individual degree of the monomials, to take the
/// elementwise maximum of two monomial, we can use the lcm function.
pub fn create_kron_data<'a, P>(poly_a: &Poly<'a, P>, poly_b: &Poly<'a, P>) -> KroneckerData<'a, P>
where
    P: PolyRing,
{
    let max_degs_a = poly_a
        .terms
        .iter()
        .fold(P::Mon::zero(), |acc, term| acc.lcm(&term.mon));
    let max_degs_b = poly_b
        .terms
        .iter()
        .fold(P::Mon::zero(), |acc, term| acc.lcm(&term.mon));
    KroneckerData {
        sub_values: (max_degs_a + max_degs_b).into_iter().collect(),
        ring: poly_a.ring.unwrap()
    }
}

/// Creates a Univariate polynomial out of a multivariate one and instructions for its
/// Kronecker substitution
///
/// FIXME We clone the monomial each time. This is very inefficient! If we use the
/// MultiIndex monomial struct, then we would easily avoid this because we could implement
/// the Iterator trait, but because of more exotic monomial representations we can't
/// implement Iterator on all of them.
pub fn to_univar<'a, P: PolyRing>(
    poly_a: &Poly<'a, P>,
    kron_data: KroneckerData<P>,
) -> PolyU<'a, P::Coeff> {
    let new_terms = poly_a
        .terms
        .iter()
        .map(|term| {
            let new_deg = term.mon.clone().into_iter().zip(kron_data.sub_values.iter()).map(|(a, b)| a * *b)
                .sum();
            Term {
                coeff: term.coeff,
                mon: UniIndex(new_deg),
            }

        })
        .collect();
    Poly::from_terms(new_terms, None)
}

/// Uses Kronecker substitution to convert the multivariate polynomial to a univariate one
pub fn to_multivar<'a, P: PolyRing>(
    poly_a: PolyU<'a, P::Coeff>,
    kron_data: KroneckerData<'a, P>,
) -> Poly<'a, P> {
    let terms = poly_a.terms.iter().map(|term| {
                let mut deg_usize = term.mon.0;
                let new_monomial = kron_data.sub_values.iter().map(|d| {
                    let r = deg_usize.rem_euclid(*d);
                    deg_usize = deg_usize.div_euclid(*d);
                    r
                    }).collect();
        Term::new(term.coeff, new_monomial) 
    }).collect();
    Poly::from_terms(terms, Some(kron_data.ring))

}

impl<'a, P: PolyRing> Poly<'a, P> {
    pub fn kronecker(
        poly_a: Self,
        poly_b: Self,
    ) -> (PolyU<'a, P::Coeff>, PolyU<'a, P::Coeff>, KroneckerData<'a, P>) {
        let kron_data = create_kron_data(&poly_a, &poly_b);
        let poly_a_uni = to_univar(&poly_a, kron_data.clone());
        let poly_b_uni = to_univar(&poly_b, kron_data.clone());

        (poly_a_uni, poly_b_uni, kron_data)
    }
}
