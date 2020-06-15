pub mod f4;
pub mod groebner_basis;
use crate::algebras::polyring::*;

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
    pub fn add(&mut self, item: Poly<'a, P>) {
        if !self.gens.contains(&item) && !item.is_zero() {
            self.gens.push(item)
        }
    }
    pub fn num_gens(&self) -> usize {
        self.gens.len()
    }
}

impl<'a, P: FPolyRing> Poly<'a, P> {
    fn reduce(&self, g: &Ideal<'a, P>) -> Option<Self> {
        // Obtain a vector of references to the polynomials 
        // as required by the euclid_div_multi function
        let ref_vec = g.gens.iter().collect();
        Some(self.euclid_div_multi(&ref_vec)?.1)
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
    gens: Vec<Term<P>>,
    original: Option<&'b Ideal<'a, P>>,
}

use crate::algebras::*;

impl<'a, 'b, P: FPolyRing> MonomialIdeal<'a, 'b, P> {
    // From a vector of terms
    pub fn new(gens: Vec<Term<P>>) -> Self {
        // Take out any zeros
        MonomialIdeal {
            gens: gens.into_iter().filter(|t| *t != Term::zero()).collect(),
            original: None,
        }
    }

    // From a polynomial idea. Keeps a reference to that ideal.
    pub fn from(poly_ideal: &'b Ideal<'a, P>) -> Self {
        MonomialIdeal {
            gens: poly_ideal.gens.iter().map(|a| a.lt().clone()).collect(),
            original: Some(poly_ideal),
        }
    }

    // Add an element to the ideal
    pub fn add(&mut self, item: Term<P>) {
        if !self.is_in(&item) && item != Term::zero() {
            self.gens.push(item)
        }
    }

    // Check if a monomial is contained in the ideal
    pub fn is_in(&self, term: &Term<P>) -> bool {
        self.gens.iter().any(|t| t.divides(term).unwrap())
    }

    // If term is in the ideal it will return the polynomial that divides it.
    // Will return a none if it isn't in the ideal or the monomial ideal wasn't
    // derived from a polynomial ideal
    pub fn get_poly(&self, term: &Term<P>) -> Option<&Poly<'a, P>> {
        match self.original {
            Some(ideal) => {
                for (lm, poly) in izip!(self.gens.iter(), ideal.gens.iter()) {
                    if lm.divides(term).unwrap() {
                        return Some(&poly);
                    }
                }
                None
            }
            None => None,
        }
    }
}
