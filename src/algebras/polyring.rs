use crate::polyu::*;
use crate::algebras::*;

impl<T: Group> AlgAdd for PolyU<T> {
    fn add(&self, other: &Self) -> Self {
        PolyU::elementwise_binary(&self, other, |a, b| a.add(&b))
    }
}

impl<T: Group> AlgAddAssign for PolyU<T> {
    fn add_ass(&mut self, other: &Self) {
        PolyU::elementwise_binary_ass(self, other, |a, b| a.add(&b))
    }
}

impl<T: Group> AlgSub for PolyU<T> {
    fn sub(&self, other: &Self) -> Self {
        PolyU::elementwise_binary(self, other, |a, b| a.sub(&b))
    }
}

impl<T: Group> AlgSubAssign for PolyU<T> {
    fn sub_ass(&mut self, other: &Self) {
        PolyU::elementwise_binary_ass(self, other, |a, b| a.sub(&b))
    }
}

impl<T: Group> AlgNeg for PolyU<T> {
    fn neg(&self) -> Self {
        let terms = self.terms.iter().map(|m| m.neg()).collect();
        PolyU::from_monomials(self.symb.clone(), terms).unwrap()
    }
}

impl<T: Group> Group for PolyU<T> {
    fn zero() -> Self {
        PolyU { symb: None, terms: vec![Monomial::zero()] }
    }
}


impl<T: Ring> AlgMul for PolyU<T> {
    fn mul(&self, other: &Self) -> Self {
        let mut result = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                result[a.deg + b.deg].add_ass(&a.coeff.mul(&b.coeff));
            }
        }

        PolyU::from_coeff( self.symb.clone(), result ).unwrap()
    }
}

impl<T: Ring> AlgMulAssign for PolyU<T>{
    fn mul_ass(&mut self, other: &Self) {
        let mut result = self.mul(other);
        self.terms.clear();
        self.terms.append(&mut result.terms)
    }
}

impl<T: Ring> Ring for PolyU<T> {
    fn one() -> Self {
        PolyU { symb: None, terms: vec![Monomial::one()] }
    }
    type BaseRing = T;

    fn scale(&self, scalar: Self::BaseRing) -> Self {
        if scalar == <T>::zero() {
            return PolyU::zero()
        }

        let result = self.terms.iter()
                         .map(|x|
                            Monomial::new(x.coeff.mul(&scalar), x.deg)
                        ).collect();

        PolyU::from_monomials(
            self.symb.clone(),
            result,
        ).unwrap()

    }
    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        // TODO need to compress if the scalar was zero
        for x in self.terms.iter_mut() {
            // Yes I'm doing mul_ass instead of scale_ass
            x.coeff.mul_ass(&scalar);
        }
    }
}

impl<T: Group> PolyU<T> {
    // TODO This whole thing needs to be redone
    fn elementwise_binary<F>(polya: &PolyU<T>, polyb: &PolyU<T>, func: F) -> PolyU<T>
    where
        F: Fn(T, T) -> T
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
                    let a = func(x.coeff.clone(), y.coeff.clone());
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

    // Pretty hacky way but okay
    fn elementwise_binary_ass<F>(polya: &mut PolyU<T>, polyb: &PolyU<T>, func: F)
    where
        F: Fn(T, T) -> T
    {
        let mut result = PolyU::elementwise_binary(polya, polyb, func);

        polya.terms.clear();
        polya.terms.append(&mut result.terms);
    }
}
