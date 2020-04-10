use crate::algebras::*;
use std::marker::PhantomData;

// <><><><><><><><><><> Poly Ring Domain <><><><><><><><><><> //
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PRDomain<T: ScalarRing, U: Variate> {
    pub(crate) vars: U,
    coeff_type: PhantomData<T>,
}


impl<T: ScalarRing, U: Variate> PRDomain<T, U> {
    pub fn new(vars: U) -> Self { PRDomain { vars, coeff_type: PhantomData } }
}

type SymbType = String;

// <><><><><><><><><><> Polynomial <><><><><><><><><><> //
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Poly<'a, T: ScalarRing, U: Variate> {
    pub(crate)    ls : T,
    pub(crate) terms : Vec<Term<T>>,
    pub(crate)  ring : &'a PRDomain<T, U>,
}

impl<'a, T: ScalarRing, U: Variate> Poly<'a, T, U> {
    // TODO Probably should make a checked version where I sort the list and remove duplicates
    pub fn from_terms_unchecked(terms: Vec<Term<T>>, ring: &'a PRDomain<T, U>) -> Poly<'a, T, U> {
        Poly { ls: <T>::one(), terms, ring }
    }
}


// <><><><><><><><><><> Term <><><><><><><><><><> //
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Term<T: ScalarRing> {
    pub coeff : T,
    pub deg   : TermIndex,
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct TermIndex ( pub(crate) u64 );

impl Group for TermIndex {
    fn zero() -> Self { TermIndex(0) }

    fn add(&self, other: &Self) -> Self {
        TermIndex( self.0 + other.0 )
    }
    fn sub(&self, other: &Self) -> Self {
        TermIndex( self.0 - other.0 )
    }
    fn neg(&self) -> Self {
        panic!("No thank you")
    }
}

impl<T: ScalarRing> Term<T> {
    
    pub fn new(coeff: T, deg: TermIndex) -> Self {
        Term { coeff, deg }
    }
    pub fn zero() -> Self { Term::new( <T>::zero(), TermIndex(0) ) }

    pub fn scale(&self, scalar: &T) -> Self {
        Term::new( self.coeff.mul(scalar), self.deg )
    } 
}

// <><><><><><><><><><> General Polynomial Functions <><><><><><><><><><> //
impl<'a, T: ScalarRing, U: Variate> Poly<'a, T, U> {

    // Only want these two arrays to be incremented at the same time
    pub fn push(&mut self, (coeff, deg): (T, TermIndex)) {
        self.terms.push(Term::new(coeff, deg));
    }

    pub fn get(&self, i: usize) -> Option<&Term<T>> {
        self.terms.get(i)
    }

    // A get_unchecked implementation
    pub(crate) fn get_uc(&self, i: usize) -> &Term<T> {
        unsafe{
            self.terms.get_unchecked(i)
        }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<T> {
        match self.num_terms() {
            0 => Term::zero(),
            n => {
                let lt = self.get(n-1).unwrap();
                lt.scale(&self.ls)
            }
        }
    }

    pub fn num_terms(&self) -> usize     { self.terms.len() }
    pub fn lc(&self)        -> T         { self.lt().coeff }
    pub fn lm(&self)        -> TermIndex { self.lt().deg }
    pub fn deg(&self)       -> usize     { self.ring.vars.tdeg(&self.lm()) }
}


pub trait Variate: Eq + PartialEq + Clone + std::fmt::Debug {
    fn cmp(&self, a: &TermIndex, b: &TermIndex) -> Ordering;
    fn tdeg(&self, index: &TermIndex) -> usize;
    fn zero() -> TermIndex;
}

impl<'a, T: ScalarRing, U: Variate> PolyRing for Poly<'a, T, U> {
    type BaseRing = T;

    fn is_zero(&self) -> bool {
        self.num_terms() == 1 && self.get(0).unwrap().coeff == <T>::zero()
    }

    fn add(&self, other: &Self) -> Self {
        Poly::elementwise_add(&self, other)
    }

    // TODO the negation here is VERY expensive because it is cloning everything 
    // The best solution I can think of is to use smart pointers, and hold things
    // like the lead scalar in the smart pointer
    fn sub(&self, other: &Self) -> Self {
        Poly::elementwise_add(self, &other.neg())
    }

    fn neg(&self) -> Self {
        let mut res = self.clone();
        res.ls = res.ls.neg();
        res
    }

    fn zero(&self) -> Self {
        Poly::from_terms_unchecked(vec![Term::zero()], self.ring)
    }

    fn is_one(&self) -> bool {
        self.num_terms() == 1 && self.get(0).unwrap().coeff == <T>::one()
    }

    fn one(&self) -> Self {
        Poly::from_terms_unchecked(vec![Term::zero()], self.ring)
    }

    fn scale(&self, scalar: Self::BaseRing) -> Self {
        if scalar == <Self::BaseRing>::zero() {
            Poly::zero(&self)
        } else {
            let mut result = self.clone();
            result.ls.mul_ass(&scalar);
            result
        }
    }

    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        if scalar == <Self::BaseRing>::zero() {
            self.clone_from(&Poly::zero(&self));
        } else {
            self.ls.mul_ass(&scalar);
        }
    }
}

use std::cmp::Ordering;

impl<'a, T: ScalarRing, U: Variate> Poly<'a, T, U> {
    
    fn elementwise_add(polya: &Self, polyb: &Self) -> Self {

        if polya.ring != polyb.ring {
            panic!("Try to add two poly with different rings I don't know how to do that yet")
        }

        // Knowing which polynomial is bigger is advantages here
        let (smol, bigg) = match polya.ring.vars.cmp(&polya.lm(), &polyb.lm()) {
            Ordering::Less => (&polya, &polyb),
            _              => (&polyb, &polya),
        };

        // Note: with_capacity allocates memory, should be deallocate some after we
        //       finished?
        let mut res = Vec::with_capacity(bigg.num_terms() + smol.num_terms());

        let mut i = 0;
        let mut j = 0;

        while i < smol.num_terms() {
            let (a, b) = (smol.get_uc(i), bigg.get_uc(j));
            res.push( match polya.ring.vars.cmp(&a.deg, &b.deg) { // Compare their degrees
                Ordering::Less    => { i += 1; a.scale(&smol.ls) },
                Ordering::Greater => { j += 1; b.scale(&bigg.ls) },
                Ordering::Equal   => {
                    i += 1; j += 1;
                    let c = a.coeff.mul(&smol.ls).add(&b.coeff.mul(&bigg.ls));
                    if c == <T>::zero() {
                        continue
                    }
                    Term { coeff: c, deg: a.deg }
                },
            })
        }

        // Append any remaining terms to the result vector
        for k in j..bigg.num_terms() {
            res.push(bigg.get_uc(k).clone())
        }

        Poly::from_terms_unchecked(res, polya.ring)
    }
}

// <><><><><><><><> Iterators <><><><><><><><> //
// pub(crate) struct PolyIter<'a, T: ScalarRing, U: Variate> {
//     data: &'a Poly<'a, T, U>,
//     index: usize,
// }

// impl<'a, T: ScalarRing, U: Variate> Iterator for PolyIter<'a, T, U> {

//     type Item = Term<T>;

//     fn next(&mut self) -> Option<Self::Item> {
//         self.data.terms.get(self.index).map(
//             |mono| {
//                 self.index += 1;
//                 *mono
//             }
//         )
//     }
// }

// impl<'a, T: ScalarRing, U: Variate> Poly<'a, T, U> {
//     pub(crate) fn iter(&'a self) -> PolyIter<'a, T, U> {
//         PolyIter {
//             data: &self,
//             index: 0
//         }
//     }
// }