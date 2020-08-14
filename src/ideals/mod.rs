use crate::algebras::*;
use crate::algebras::polyring::*;
use num_traits::Zero;

pub mod f4;
pub mod groebner_basis;

#[derive(Debug, Clone)]
pub struct Ideal<'a, P: PolyRing> {
    pub gens: Vec<Poly<'a, P>>,
}

impl<'a, P: PolyRing> Ideal<'a, P> {
    pub fn new(gens: Vec<Poly<'a, P>>) -> Self {
        Ideal { 
            gens: gens.into_iter().filter(|t| !t.is_zero()).collect()
        }
    }
    /// Returns a reference to the item if it was successfully added
    /// Returns None if not
    pub fn add(&mut self, item: Poly<'a, P>) -> Option<&Poly<'a, P>> {
        if !self.gens.contains(&item) && !item.is_zero() {
            self.gens.push(item);
            self.gens.last()
        } else {
            None
        }
    }
    pub fn num_gens(&self) -> usize {
        self.gens.len()
    }
}

use std::fmt;

impl<'a, P: PolyRing> fmt::Display for Ideal<'a, P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Ideal (\n {} )", 
            self.gens.iter()
                .map(|p| format!("{}", p))
                .collect::<Vec<String>>()
                .join("\n")
        )
    }
}


// Want to take ownership so that the ideal doesn't change.
// If you want the ideal again then you need to get it from the destructor
#[derive(Debug, Clone)]
pub struct MonomialIdeal<'a, 'b, P: FPolyRing> {
    gens: Vec<P::Mon>,
    original: Option<&'b Ideal<'a, P>>,
}

// Should return the minimum generating set of monomials
fn simplify<P: PolyRing>(input: Vec<P::Mon>) -> Vec<P::Mon> {

    unimplemented!()
    // let n = input.len();
    // let res = Vec::with_capacity(n);
    // for i in 0..n {
    //     let is_minimal = 0..n.filter(|k| k != i).all(|j| !input[i].divides(input[j]).unwrap());

    //     if is_minimal {
    //         res.push(input[i])
    //     }
    // }
}

#[cfg(test)]
mod test {
    #[test]
    fn is_simplified() {

        // let mon_ideal = MonomialIdeal::new(Monomial

    }
}


impl<'a, 'b, P: FPolyRing> MonomialIdeal<'a, 'b, P> {
    pub fn new(gens: Vec<P::Mon>) -> Self {
        // Take out any zeros
        MonomialIdeal {
            gens: gens.into_iter().filter(|t| !t.is_zero()).collect(),
            original: None,
        }
    }
    /// Has the problem that it might include repeated elements if more than one polynomial
    /// has the same leading monomial. Might have to resolve that some time
    pub fn from(poly_ideal: &'b Ideal<'a, P>) -> Self {
        MonomialIdeal {
            gens: poly_ideal.gens.iter().map(|a| a.lm().clone()).collect(),
            original: Some(poly_ideal),
        }
    }
    pub fn add(&mut self, item: P::Mon) {
        if !self.is_in(&item) && item != P::Mon::zero() {
            self.gens.push(item)
        }
    }

    pub fn is_in(&self, term: &P::Mon) -> bool {
        self.gens.iter().any(|t| term.divides(t).unwrap())
    }

    /// If term is in the ideal it will return the element that corresponds to it
    /// This could be done faster
    pub fn get_poly(&self, term: &P::Mon) -> Option<&Poly<'a, P>> {
        match self.original {
            Some(ideal) => {
                for (lm, poly) in izip!(self.gens.iter(), ideal.gens.iter()) {
                    // println!("lm = {:?}", lm);
                    // println!("term = {:?}", term);
                    // println!("poly = {}", poly);
                    if term.divides(lm).unwrap() {
                        return Some(&poly);
                    }
                }
                None
            },
            None => None
        }
    }
}
