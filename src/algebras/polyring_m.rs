use crate::polyu::*;
use crate::error::*;
use crate::algebras::*;

impl<T: PolyRing> Ring for PolyM<T> {
    type BaseRing = <T>::BaseRing;
    fn is_poly() -> bool { true }

    fn zero() -> Self {
        PolyM::from_coeff( None, vec![<T>::zero()]).unwrap()
    }
    fn add(&self, other: &Self) -> Self {
        PolyM::elementwise_binary(&self, other, |a, b| a.add(&b))
    }

    fn sub(&self, other: &Self) -> Self {
        PolyM::elementwise_binary(&self, other, |a, b| a.sub(&b))
    }
    fn neg(&self) -> Self {
        let mut result = self.clone();
        result.lead_scalar = result.lead_scalar.neg();
        result
    }
    fn mul(&self, other: &Self) -> Self {
        let mut result = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                result[a.deg + b.deg] = result[a.deg + b.deg].add(&a.coeff.mul(&b.coeff));
            }
        }

        PolyM::from_coeff( self.symb.clone(), result ).unwrap()
    }

    fn one() -> Self {
        PolyM::from_coeff( None, vec![<T>::one()]).unwrap()
    }
}

use std::fmt::Debug;

impl<T: PolyRing + Debug> PolyRing for PolyM<T> {

    fn set_symb(&mut self, new_symb: SymbType) {
        self.symb = new_symb;
    }

    fn get_symb(&self) -> SymbType {
        self.symb
    }

    fn scale(&self, scalar: Self::BaseRing) -> Self {

        if scalar == <Self::BaseRing>::zero() {
            return Poly::zero()
        }

        let mut result = self.clone();
        result.lead_scalar.mul_ass(&scalar);
        result
    }

    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        if scalar == <Self::BaseRing>::zero() {
            self.terms.clear()
        } else {
            self.lead_scalar.mul_ass(&scalar);
        }
    }
}

impl<T: PolyRing> PolyM<T> {

    pub fn from_coeff(symb: SymbType, coeffs: Vec<T>) -> Result<PolyM<T>, PolyErr> {
        // Converts into a PolyU type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument

        let mut terms = Vec::new();
        for (i, c) in coeffs.into_iter().enumerate() {
            if c != <T>::zero() {
                terms.push(MonomialM::<T>::new(c, i));
            }
        }

        Ok(Poly::from_monomials(symb, terms).unwrap())
    }

    // TODO This is painful to read
    fn elementwise_binary<F>(polya: &PolyM<T>, polyb: &PolyM<T>, func: F) -> PolyM<T>
    where
        F: Fn(T, T) -> T
    {
        if polya.terms.is_empty() { 
            return polyb.clone()
        } else if polyb.terms.is_empty() {
            return polya.clone()
        }
        // From now one we assume that the vectors are nonempty
        if polya.symb != polyb.symb {
            let coerce = PolyM::from_coeff(polya.symb, vec![polyb]);
            let first = polya.terms[0].coeff.add(&polyb.terms);
        }
        let (smol, bigg) =
            if polya.deg() > polyb.deg() {
                (&polyb.terms, &polya.terms)
            } else {
                (&polya.terms, &polyb.terms)
            };

        let mut result: Vec<MonomialM<T>> = Vec::with_capacity(bigg.len());

        let mut i = 0;
        let mut j = 0;

        while i < smol.len() {
            match (smol[i].clone(), bigg[j].clone()) {
                (x, y) if x.deg <  y.deg => {result.push(x.clone()); i += 1},
                (x, y) if x.deg >  y.deg => {result.push(y.clone()); j += 1}, 
                (x, y) if x.deg == y.deg => {
                    let a = func(x.coeff.clone(), y.coeff.clone());
                    if a != <T>::zero() {
                        result.push(Monomial::new(a, x.deg));
                    }
                    i += 1; j += 1;
                },
                _ => unreachable!(),
            };
        }

        // Append any remaining terms to the result vector
        for j in i..bigg.len() {
            result.push(bigg[j].clone())
        }

        if result.is_empty() {
            result.push(Monomial::zero())
        }

        Poly::from_monomials(polya.symb, result).unwrap()
    }
}