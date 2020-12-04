use crate::algebras::polyring::*;
use crate::algebras::*;
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
            gens: gens.into_iter().filter(|t| !t.is_zero()).collect(),
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
        write!(
            f,
            "Ideal (\n {} )",
            self.gens
                .iter()
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
                    if term.divides(lm).unwrap() {
                        return Some(&poly);
                    }
                }
                None
            }
            None => None,
        }
    }
}

// #[cfg(test)]
// mod test {

//     use crate::algebras::real::*;
//     use crate::parse::MyFromStr;
//     use crate::polym::*;
//     use generic_array::typenum::{U2, U3};

//     use crate::algebras::polyring::*;
//     use crate::display::*;
//     use crate::ideals::*;
//     use num_traits::Zero;
//     #[test]
//     fn is_in() {
//         let ideal = MonomialIdeal::new((0..100).map(|i| {
//             Term::new(ZZ(1), 
//                       MultiIndex::new((0..U3::to_usize()).map(|i| rand()).collect::<GenericArray<usize, U3>>()))
//         }).collect();
//     }
// }
