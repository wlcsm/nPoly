use crate::polyu::*;
use crate::algebras::*;

impl<T: ScalarRing + std::fmt::Debug> Ring for PolyU<T> {

    type BaseRing = T;

    fn is_poly() -> bool { true }

    fn add(&self, other: &Self) -> Self {
        PolyU::elementwise_add(&self, other)
    }

    // TODO the negation here is VERY expensive because it is cloning everything 
    // The best solution I can think of is to use smart pointers, and hold things
    // like the lead scalar in the smart pointer
    fn sub(&self, other: &Self) -> Self {
        PolyU::elementwise_add(self, &other.neg())
    }

    fn neg(&self) -> Self {
        let mut res = self.clone();
        res.terms.lead_scalar = res.terms.lead_scalar.neg();
        res
    }

    fn zero() -> Self {
        PolyU::from_coeff(None, vec![<T>::zero()]).unwrap()
    }

    /// Standard O(n^2) multiplication
    fn mul(&self, other: &Self) -> Self {

        let mut acc = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for (coeff_a, deg_a) in self.terms.iter() {
            for (coeff_b, deg_b) in other.terms.iter() {
                acc[deg_a + deg_b].add_ass(&coeff_a.mul(&coeff_b));
            }
        }

        PolyU::from_coeff(self.symb.clone(), acc).unwrap()
    }

    fn one() -> Self {
        PolyU::from_coeff(None, vec![<T>::one()]).unwrap()
    }
}

impl<T: ScalarRing> PolyRing for PolyU<T> {

    // These functions clear all the information if the scalar is zero.
    // Thats good but it can take a little bit of time perhaps?

    fn scale(&self, scalar: Self::BaseRing) -> Self {
        if scalar == <Self::BaseRing>::zero() {
            PolyU::zero()
        } else {
            let mut result = self.clone();
            result.terms.lead_scalar.mul_ass(&scalar);
            result
        }
    }

    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        if scalar == <Self::BaseRing>::zero() {
            self.clone_from(&PolyU::zero());
        } else {
            self.terms.lead_scalar.mul_ass(&scalar);
        }
    }
}

use std::cmp::Ordering;

impl<T: ScalarRing> PolyU<T> {
    
    fn elementwise_add(polya: &PolyU<T>, polyb: &PolyU<T>) -> PolyU<T> {

        let (smol, bigg) = match polya.deg().cmp(&polyb.deg()) {
            Ordering::Less => (&polya.terms, &polyb.terms),
            _              => (&polyb.terms, &polya.terms),
        };

        let mut res = Monomials::with_capacity(bigg.len() + smol.len());

        let mut i = 0;
        let mut j = 0;

        while i < smol.len() {
            let (a, b) = (smol.get_uc(i), bigg.get_uc(j));
            res.push( match a.1.cmp(&b.1) { // Compare their degrees
                Ordering::Less    => { i += 1; (a.0.mul(&smol.lead_scalar), a.1) },
                Ordering::Greater => { j += 1; (b.0.mul(&bigg.lead_scalar), b.1) },
                Ordering::Equal   => {
                    i += 1; j += 1;
                    let c = a.0.mul(&smol.lead_scalar).add(&b.0.mul(&bigg.lead_scalar));
                    if c == <T>::zero() {
                        continue
                    }
                    (c, a.1)
                },
            })
        }

        // Append any remaining terms to the result vector
        for k in j..bigg.len() {
            res.push(bigg.get_uc(k))
        }

        PolyU::from_terms(polya.symb.clone(), res).unwrap()
    }
}
