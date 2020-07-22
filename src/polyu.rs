use crate::algebras::polyring::*;
use crate::algebras::*;
use num_traits::Zero;
use std::cmp::Ordering;

// Shorthand
pub type PolyU<'a, R> = Poly<'a, PRDomain<R, UniVarOrder>>;

pub trait PolyRingUni: PolyRing<Ord=UniVarOrder, Mon=UniIndex> {}

impl<R: ScalarRing> PolyRingUni for PRDomain<R, UniVarOrder> {}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct UniVarOrder();

impl MonOrd for UniVarOrder {
    type Index = UniIndex;

    fn cmp(a: &UniIndex, b: &UniIndex) -> Ordering {
        a.0.cmp(&b.0)
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

use std::ops;

// Since UniIndex implements Copy, a default implementation for MyAddMonoid
// is provided
impl ClosedAdd for UniIndex {}
impl_op_ex!(+ |a: &UniIndex, b: &UniIndex| -> UniIndex { UniIndex(a.0 + b.0) });
impl_op_ex!(+= |a: &mut UniIndex, b: &UniIndex| { a.0 += b.0 });

// impl MyAddMonoid for UniIndex {
//     fn ref_add(&self, other: &Self) -> Self {
//         UniIndex(self.0 + other.0)
//     }
// }

use std::cmp::{max, min};
use typenum::U1;
impl VarNumber for U1 {}

impl Monomial for UniIndex {

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

    fn div(&self, other: &Self) -> Option<Self> {
        if self.0 < other.0 {
            None
        } else {
            Some(UniIndex(self.0 - other.0))
        }
    }

    fn lex(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }

    // fn sub(&self, other: &Self) -> Option<Self> {
    //     if other.0 <= self.0 {
    //         Some(UniIndex(self.0 - other.0))
    //     } else {
    //         None
    //     }
    // }

    // fn divides(&self, other: &Self) -> Option<bool> {
    //     Some(self.0 <= other.0)
    // }
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
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| Term::new(c, UniIndex(i)))
            .collect();

        Poly::from_terms_unchecked(terms, Some(ring))
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
                // I still don't know why I can't use method syntax for the tdeg calls...
                acc[d_a.0 + <UniIndex>::tot_deg(d_b)] += *c_a * *c_b;
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
