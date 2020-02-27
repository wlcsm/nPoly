use crate::polyu::*;
use crate::algebras::*;

impl<T: ScalarRing + std::fmt::Debug> Ring for PolyU<T> {
    type BaseRing = T;
    fn is_poly() -> bool { true }

    fn add(&self, other: &Self) -> Self {
        PolyU::elementwise_add(&self, other)
    }
    // TODO the negation here is VERY expensive because it is cloning everything 
    fn sub(&self, other: &Self) -> Self {
        PolyU::elementwise_add(self, &other.neg())
    }
    fn neg(&self) -> Self {
        let terms = self.terms.iter().map(|m| m.neg()).collect();
        PolyU::from_monomials(self.symb.clone(), terms).unwrap()
    }
    fn zero() -> Self {
        PolyU { symb: None, lead_scalar: <T>::one(), terms: vec![Monomial::zero()] }
    }
    fn mul(&self, other: &Self) -> Self {
        let mut acc = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                acc[a.deg + b.deg].add_ass(&a.coeff.mul(&b.coeff));
            }
        }
        // Convert everything to monomials
        let mut result = Vec::with_capacity(acc.len());
        for (i, el) in acc.into_iter().enumerate() {
            if el != <T>::zero() {
                result.push(Monomial::new(el, i));
            }
        }

        PolyU::<T> { symb: self.symb.clone(), lead_scalar: <T>::one(), terms: result }
    }

    fn one() -> Self {
        PolyU { symb: None, lead_scalar: <T>::one(), terms: vec![Monomial::one()] }
    }
}

impl<T: ScalarRing> PolyRing for PolyU<T> {
    fn scale(&self, scalar: Self::BaseRing) -> Self {
        if scalar == <Self::BaseRing>::zero() {
            PolyU::zero()
        } else {
            let mut result = self.clone();
            result.lead_scalar.mul_ass(&scalar);
            result
        }
    }
    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        if scalar == <Self::BaseRing>::zero() {
            self.terms.clear();
            self.lead_scalar = <Self::BaseRing>::one();
        } else {
            self.lead_scalar.mul_ass(&scalar);
        }
    }
}

impl<T: ScalarRing> PolyU<T> {
    // TODO This whole thing needs to be redone
    fn elementwise_add(polya: &PolyU<T>, polyb: &PolyU<T>) -> PolyU<T>
    {
        let (smol, bigg) =
            if polya.deg() > polyb.deg() {
                (&polyb.terms, &polya.terms)
            } else {
                (&polya.terms, &polyb.terms)
            };

        let mut result: Vec<Monomial<T>> = Vec::with_capacity(bigg.len());

        let mut i = 0;
        let mut j = 0;

        while i < smol.len() {
            match (smol[i].clone(), bigg[j].clone()) {
                (x, y) if x.deg <  y.deg => {result.push(x.clone()); i += 1},
                (x, y) if x.deg >  y.deg => {result.push(y.clone()); j += 1}, 
                (x, y) if x.deg == y.deg => {
                    i += 1; j += 1;
                    let a = x.coeff.add(&y.coeff);
                    if a == <T>::zero() {
                        result.push(Monomial::new(a, x.deg));
                    }
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

        PolyU::from_monomials(polya.symb.clone(), result).unwrap()
    }
}
