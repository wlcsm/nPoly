use crate::algebras::*;
use std::cmp::Ordering;
use std::marker::PhantomData;
use generic_array::*;
use std::fmt::Debug;

// <><><><><><><><><><> Poly Ring Domain <><><><><><><><><><> //
pub(crate) trait PolyRing: Eq + PartialEq + Clone {
    type Coeff: ScalarRing;
    type Var: Variate;
    fn symb(&self) -> Vec<String>;
}

// The reason I have kept this generic for now (when at the moment there isn't any immediate
// need to) is because I plan on implementing quotient rings for which I will need access to generic parameters
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PRDomain<T: ScalarRing, U: Variate> {
    pub(crate) var: GenericArray<usize, U::NumVar>,
    coeff_type: PhantomData<T>,
    index_type: PhantomData<U>,
}

impl<T: ScalarRing, U: Variate> PolyRing for PRDomain<T, U> {
    type Coeff = T;
    type Var = U;
    fn symb(&self) -> GenericArray<usize, U::NumVar> { self.var }
}

impl<T: ScalarRing, U: Variate> PRDomain<T, U> {
    pub fn new(vars: GenericArray<usize, U::NumVar>) -> Self { 
        PRDomain { 
            var: vars,
            coeff_type: PhantomData,
            index_type: PhantomData,
        }
    }
}

// <><><><><><><><><><> Polynomial <><><><><><><><><><> //
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Poly<'a, P: PolyRing> {
    pub(crate)    ls : P::Coeff,
    pub(crate) terms : Vec<Term<P>>,
    pub(crate)  ring : &'a P,
}

impl<'a, P: PolyRing> Poly<'a, P> {
    // TODO Probably should make a checked version where I sort the list 
    // and remove duplicates
    pub fn from_terms_unchecked(terms: Vec<Term<P>>, ring: &'a P) -> Poly<'a, P> {
        Poly { ls: <P::Coeff>::one(), terms, ring }
    }
}

// <><><><><><><><><><> Term <><><><><><><><><><> //
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Term<P: PolyRing> {
    pub coeff : P::Coeff,
    pub deg   : <P::Var as Variate>::Index,
}


impl<P: PolyRing> Term<P> {
    
    pub fn new(coeff: P::Coeff, deg: <P::Var as Variate>::Index) -> Self {
        Term { coeff, deg }
    }
    pub fn zero() -> Self { Term::new( <P::Coeff>::zero(), <P::Var as Variate>::Index::zero() ) }

    pub fn scale(&self, scalar: &P::Coeff) -> Self {
        Term::new( self.coeff.mul(scalar), self.deg )
    } 
}

// <><><><><><><><><><> General Polynomial Functions <><><><><><><><><><> //
impl<'a, P: PolyRing> Poly<'a, P> {

    // Only want these two arrays to be incremented at the same time
    pub fn push(&mut self, (coeff, deg): (P::Coeff, <P::Var as Variate>::Index)) {
        self.terms.push(Term::new(coeff, deg));
    }

    pub fn get(&self, i: usize) -> Option<&Term<P>> {
        self.terms.get(i)
    }

    // A get_unchecked implementation
    pub(crate) fn get_uc(&self, i: usize) -> &Term<P> {
        unsafe{
            self.terms.get_unchecked(i)
        }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<P> {
        match self.num_terms() {
            0 => Term::zero(),
            n => {
                let lt = self.get(n-1).unwrap();
                lt.scale(&self.ls)
            }
        }
    }

    pub fn num_terms(&self) -> usize     { self.terms.len() }
    pub fn lc(&self)        -> P::Coeff  { self.lt().coeff }
    pub fn deg(&self)       -> usize     { <P::Var>::tdeg(&self.lm()) }
    pub fn lm(&self)        -> <P::Var as Variate>::Index    { self.lt().deg }
}

pub trait IndexTrait: Zero + Clone + Eq + Debug {
    fn deg(&self) -> usize;
    fn divides(&self, other: Self) -> bool;
}

pub trait Variate: Clone + Eq + Debug {
    type NumVar: ArrayLength<usize>;
    type Index: IndexTrait;

    fn cmp(a: &Self::Index, b: &Self::Index) -> Ordering;
    fn tdeg(index: &Self::Index) -> usize;
    fn zero() -> Self::Index;
}

// Basic Arithmetic operations for polynomial
impl<'a, P: PolyRing> Poly<'a, P> {

    pub fn is_zero(&self) -> bool {
        self.num_terms() == 1 && self.get(0).unwrap().coeff == <P::Coeff>::zero()
    }

    pub fn add(&self, other: &Self) -> Self {
        Poly::elementwise_add(&self, other)
    }

    // TODO the negation here is VERY expensive because it is cloning everything 
    // The best solution I can think of is to use smart pointers, and hold things
    // like the lead scalar in the smart pointer
    pub fn sub(&self, other: &Self) -> Self {
        Poly::elementwise_add(self, &other.neg())
    }

    pub fn neg(&self) -> Self {
        let mut res = self.clone();
        res.ls = res.ls.neg();
        res
    }

    pub fn zero(&self) -> Self {
        Poly::from_terms_unchecked(vec![Term::zero()], self.ring)
    }

    pub fn is_one(&self) -> bool {
        self.num_terms() == 1 && self.get(0).unwrap().coeff == <P::Coeff>::one()
    }

    pub fn one(&self) -> Self {
        Poly::from_terms_unchecked(vec![Term::zero()], self.ring)
    }

    pub fn scale(&self, scalar: P::Coeff) -> Self {
        if scalar == <P::Coeff>::zero() {
            Poly::zero(&self)
        } else {
            let mut result = self.clone();
            result.ls.mul_ass(&scalar);
            result
        }
    }

    pub fn scale_ass(&mut self, scalar: P::Coeff) {
        if scalar == <P::Coeff>::zero() {
            self.clone_from(&Poly::zero(&self));
        } else {
            self.ls.mul_ass(&scalar);
        }
    }

    pub fn elementwise_add(polya: &Self, polyb: &Self) -> Self {

        if polya.ring != polyb.ring {
            panic!("Try to add two poly with different rings I don't know how to do that yet")
        }

        // Knowing which polynomial is bigger is advantages here
        let (smol, bigg) = match <P::Var>::cmp(&polya.lm(), &polyb.lm()) {
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
            res.push( match <P::Var>::cmp(&a.deg, &b.deg) { // Compare their degrees
                Ordering::Less    => { i += 1; a.scale(&smol.ls) },
                Ordering::Greater => { j += 1; b.scale(&bigg.ls) },
                Ordering::Equal   => {
                    i += 1; j += 1;
                    let c = a.coeff.mul(&smol.ls).add(&b.coeff.mul(&bigg.ls));
                    if c == <P::Coeff>::zero() {
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