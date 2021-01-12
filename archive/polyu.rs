use crate::algebras::polyring::*;
use crate::algebras::*;
use num_traits::Zero;

// Shorthand
pub type PolyU<'a, R> = Poly<'a, PRDomain<R, 1>>;
//pub trait PolyRingUni: PolyRing<VARS = 1> {}

//impl<R: ScalarRing> PolyRingUni for PRDomain<R, 1> {}

impl std::iter::FromIterator<usize> for UniIndex {
    fn from_iter<I: IntoIterator<Item = usize>>(iter: I) -> Self {
        UniIndex(iter.into_iter().next().unwrap())
    }
}

#[derive(Debug, Eq, PartialEq, Clone, Copy, Ord, PartialOrd, Hash)]
pub struct UniIndex(pub(crate) usize);

impl Zero for UniIndex {
    fn zero() -> Self {
        UniIndex(0)
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

use std::ops::{Add, AddAssign};

impl Add for UniIndex {
    type Output = UniIndex;

    fn add(self, other: Self) -> UniIndex {
        UniIndex(self.0 + other.0)
    }
}

impl AddAssign for UniIndex {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0
    }
}

use std::cmp::{min, max};

impl EuclideanDomain for UniIndex {
    fn gcd(&self, other: &Self) -> Self { UniIndex(min(self.0, other.0)) }
    fn lcm(&self, other: &Self) -> Self { UniIndex(max(self.0, other.0)) }
    /// In the form: self / other = Some(Quotient, Remainder)
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)> {
        if self.0 < other.0 {
            Some((UniIndex(0), UniIndex(0)))
        } else {
            Some((UniIndex(self.0 - other.0), UniIndex(0)))
        }
    }
}

//impl Monomial for UniIndex {
//    type NumVar = U1;
//
//    // TODO These get and set are not good and I should look for a way around this
//    fn get(&self, _ind: usize) -> Option<&usize> {
//        Some(&self.0)
//    }
//    fn tot_deg(self: &Self) -> usize {
//        self.0
//    }
//}
//
// <><><><><><><><> Constructors <><><><><><><><> //
impl<'a, P: PolyRing> Poly<'a, P> {
    pub(crate) fn from_coeff(ring: &'a P, coeffs: Vec<R>) -> Poly<'a, P> {
        // Automatically compress the terms argument
        unimplemented!()
//        let terms = coeffs
//            .into_iter()
//            .enumerate()
//            .filter(|(_, c)| !c.is_zero())
//            .map(|(i, c)| Term::new(c, UniIndex(i)))
//            .collect();
//
//        Poly::from_terms_unchecked(terms, Some(ring))
    }
}

impl<'a, F: ScalarField> PolyU<'a, F> {
    /// Standard O(n^2) multiplication
    pub fn mul(&self, other: &Self) -> Self {
        let mut acc = vec![<F>::zero(); self.deg() + other.deg() + 1];

        for Term {
            coeff: c_a,
            mon: d_a,
        } in self.terms.iter()
        {
            for Term {
                coeff: c_b,
                mon: d_b,
            } in other.terms.iter()
            {
                unimplemented!()
                // I still don't know why I can't use method syntax for the tdeg calls...
//                acc[d_a.0 + <UniIndex>::tot_deg(d_b)] += *c_a * *c_b;
            }
        }

        Poly::from_coeff(self.ring.unwrap(), acc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::integers::ZZ;
    use crate::parse::*;

    #[test]
    fn display_test() {
        let ring = PRDomain::<ZZ, UniVarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3x^2 + 5x^98").unwrap();
        println!("{}", a);
    }
}
