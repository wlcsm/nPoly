use crate::algebras::polyring::*;
use crate::algebras::*;
use std::cmp::Ordering;

// Shorthand
pub type PolyU<'a, R> = Poly<'a, PRDomain<R, UniIndex, UnivarOrder>>;

pub trait PolyRingUni: PolyRing<Var = UniIndex, Ord = UnivarOrder> {}

impl<R: ScalarRing> PolyRingUni for PRDomain<R, UniIndex, UnivarOrder> {}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct UnivarOrder();

impl MonomialOrdering<UniIndex> for UnivarOrder {
    fn cmp(a: &UniIndex, b: &UniIndex) -> Ordering {
        a.0.cmp(&b.0)
    }
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct UniIndex(pub(crate) usize);

impl Zero for UniIndex {
    fn zero() -> Self {
        UniIndex(0)
    }
}
use std::cmp::{max, min};
use typenum::U1;
impl VarNumber for U1 {}

impl Variate for UniIndex {
    type NumVar = U1;

    // TODO These get and set are not good and I should look for a way around this
    fn get(&self, _ind: usize) -> Option<&usize> {
        Some(&self.0)
    }
    fn set(&mut self, _ind: usize, val: usize) -> Option<()> {
        self.0 = val;
        Some(())
    }
    fn tot_deg(self: &Self) -> usize {
        self.0
    }

    fn add(&self, other: &Self) -> Self {
        UniIndex(self.0 + other.0)
    }

    fn sub(&self, other: &Self) -> Option<Self> {
        if other.0 <= self.0 {
            Some(UniIndex(self.0 - other.0))
        } else {
            None
        }
    }

    fn divides(&self, other: &Self) -> Option<bool> {
        Some(self.0 <= other.0)
    }
    fn gcd(&self, other: &Self) -> Self {
        UniIndex(min(self.0, other.0))
    }
    fn lcm(&self, other: &Self) -> Self {
        UniIndex(max(self.0, other.0))
    }
}

// <><><><><><><><> Constructors <><><><><><><><> //
impl<'a, P: PolyRingUni> Poly<'a, P> {
    pub(crate) fn from_coeff(ring: &'a P, coeffs: Vec<P::Coeff>) -> Poly<'a, P> {
        // Automatically compress the terms argument
        let terms = coeffs
            .into_iter()
            .enumerate()
            .filter(|(_, c)| *c != <P::Coeff>::zero())
            .map(|(i, c)| Term::new(c, UniIndex(i)))
            .collect();

        Poly::from_terms_unchecked(terms, ring)
    }
}

impl<'a, F: Field> PolyU<'a, F> {
    /// Standard O(n^2) multiplication
    pub fn mul(&self, other: &Self) -> Self {
        let mut acc = vec![<F>::zero(); self.deg() + other.deg() + 1];

        for Term {
            coeff: c_a,
            deg: d_a,
        } in self.terms.iter()
        {
            for Term {
                coeff: c_b,
                deg: d_b,
            } in other.terms.iter()
            {
                // I still don't know why I can't use method syntax for the tdeg calls...
                acc[<UniIndex>::tot_deg(d_a) + <UniIndex>::tot_deg(d_b)].add_ass(&c_a.mul(&c_b));
            }
        }

        Poly::from_coeff(self.ring, acc)
    }
}

use std::fmt;

// Problem is that it's hard to put an ordering on the coefficients because in finite fields
// thats quite ambiguious. I need it in the "if x < 0" line
// This will eventually have to be overcome some time.
impl<'a, P: PolyRing> fmt::Display for Poly<'a, P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Because we don't want a potential "+" out the front of the first term
        if self.is_zero() {
            write!(f, "{}", <P::Coeff>::zero())
        } else {
            let mut acc: String = self.terms[0].to_str(&self.ring);

            self.terms
                .iter()
                .skip(1)
                .for_each(|x| acc.push_str(&format!(" + {}", x.to_str(&self.ring))));

            write!(f, "{}", acc)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::integers::ZZ;
    use crate::parse::*;

    #[test]
    fn display_test() {
        let ring = PRDomain::<ZZ, UniIndex, UnivarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3x^2 + 5x^98").unwrap();
        println!("{}", a);
    }
}
