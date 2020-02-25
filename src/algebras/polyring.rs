use crate::polyu::*;
use crate::error::*;
use crate::algebras::*;
use std::fmt::Debug;

impl<T: ScalarRing + Debug> Ring for Poly<T> {

    type BaseRing = T;
    fn is_poly() -> bool { true }

    fn zero() -> Self {
        Poly::from_coeff( None, vec![<T>::zero()]).unwrap()
    }

    fn add(&self, other: &Self) -> Self {
        Poly::elementwise_binary(&self, other, |a, b| a.add(&b))
    }

    fn sub(&self, other: &Self) -> Self {
        Poly::elementwise_binary(&self, other, |a, b| a.sub(&b))
    }

    fn neg(&self) -> Self {
        let result = self.clone();
        result.lead_scalar.mul_ass(&<T>::one().neg());
        result
    }

    fn mul(&self, other: &Self) -> Self {
        let mut result = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for a in self.s_terms.iter() {
            for b in other.s_terms.iter() {
                result[a.deg + b.deg].add_ass(&b.coeff.mul(&a.coeff));
            }
            for b in other.p_terms.iter() {
                result[a.deg + b.deg].add_ass(&b.coeff.scale(&<T>::one()));
            }
        }

        Poly::from_coeff( self.symb.clone(), result ).unwrap()
    }

    fn one() -> Self {
        Poly::from_coeff( None, vec![<T>::zero()]).unwrap()
    }
}

impl<T: ScalarRing + Debug> PolyRing for Poly<T> {

    fn set_symb(&mut self, new_symb: SymbType) {
        self.symb = new_symb;
    }

    fn get_symb(&self) -> SymbType {
        self.symb
    }

    fn scale(&self, scalar: T) -> Self {

        if scalar == <T>::zero() {
            Poly::zero()
        } else {
            let mut result = self.clone();
            result.lead_scalar.mul_ass(&scalar);
            result
        }
    }

    fn scale_ass(&mut self, scalar: T) {
        if scalar == <Self::BaseRing>::zero() {
            self.s_terms.clear();
            self.p_terms.clear();
            self.lead_scalar = <T>::one();
        } else {
            self.lead_scalar.mul_ass(&scalar);
        }
    }
}

impl<T: ScalarRing> Poly<T> {

    pub fn from_coeff(symb: SymbType, coeffs: Vec<T>) -> Result<Poly<T>, PolyErr> {
        // Converts into a Poly type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument

        let mut terms = Vec::new();
        for (i, c) in coeffs.into_iter().enumerate() {
            if c != <T>::zero() {
                terms.push(MonomialU::<T>::new(c, i));
            }
        }

        Ok(Poly::from_monomials(symb, terms).unwrap())
    }

    // TODO This is painful to read
    fn elementwise_binary<F>(polya: &Poly<T>, polyb: &Poly<T>, func: F) -> Poly<T>
    where
        F: Fn(T, T) -> T
    {
        let (smol, bigg) =
            if polya.deg() > polyb.deg() {
                (&polyb.terms, &polya.terms)
            } else {
                (&polya.terms, &polyb.terms)
            };

        let mut result: Vec<MonomialU<T>> = Vec::with_capacity(bigg.len());

        let mut i = 0;
        let mut j = 0;

        while i < smol.len() {
            match (smol[i].deg, bigg[j].deg) {
                (x, y) if x <  y => {result.push(smol[i].clone()); i += 1},
                (x, y) if x >  y => {result.push(bigg[j].clone()); j += 1}, 
                (x, y) if x == y => {
                    let a = func(smol[i].coeff, bigg[j].coeff);
                    if a != <T>::zero() {
                        result.push(Monomial::new(a, x));
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
