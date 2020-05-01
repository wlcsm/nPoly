use crate::algebras::*;
use std::cmp::Ordering;
use std::marker::PhantomData;
use generic_array::*;
use std::fmt::Debug;

// <><><><><><><><><><> Poly Ring Domain <><><><><><><><><><> //
pub trait PolyRing: Eq + PartialEq + Clone {
    type Coeff: ScalarRing;
    type Var: Variate;
    type Ord: MonomialOrdering<Self::Var>;
    fn symb(&self) -> Vec<usize>;
}
pub trait VarNumber: ArrayLength<usize> + Eq + PartialEq + Clone + Debug {}

// The reason I have kept this generic for now (when at the moment there isn't any immediate
// need to) is because I plan on implementing quotient rings for which I will need access to generic parameters
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PRDomain<R: ScalarRing, U: Variate, M: MonomialOrdering<U>> {
    pub(crate) var: GenericArray<usize, U::NumVar>,
    ring_parameters: PhantomData<(R, U, M)>,
}

impl<R, U, M> PolyRing for PRDomain<R, U, M>
    where R: ScalarRing, U: Variate, M: MonomialOrdering<U> {
    type Coeff = R;
    type Var = U;
    type Ord = M;
    fn symb(&self) -> Vec<usize> {
        self.var.clone().into_iter().collect()
    }
}

impl<R, U, M> PRDomain<R, U, M>
    where R: ScalarRing, U: Variate, M: MonomialOrdering<U> {
    pub fn new(vars: GenericArray<usize, U::NumVar>) -> Self { 
        PRDomain { 
            var: vars,
            ring_parameters: PhantomData,
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
    pub deg   : P::Var,
}

impl<P: PolyRing> Term<P> {
    pub fn new(coeff: P::Coeff, deg: P::Var) -> Self {
        Term { coeff, deg }
    }

    pub fn scale(&self, scalar: &P::Coeff) -> Self {
        Term::new( self.coeff.mul(scalar), self.deg.clone() )
    } 
    pub fn cmpdeg(&self, other: Self) -> Ordering {
        <P::Ord>::cmp(&self.deg, &other.deg)
    }
}
impl<P: PolyRing> Zero for Term<P> {
    fn zero() -> Self { Term::new(<P::Coeff>::zero(), <P::Var>::zero())}
}
impl<P: PolyRing> One for Term<P> {
    fn one() -> Self { Term::new(<P::Coeff>::one(), <P::Var>::zero())}
}

impl<P: PolyRing> Ring for Term<P> {

    type BaseRing = P::Coeff; 
    // Group operations
    fn add(&self, other: &Self) -> Self {
        unimplemented!()
    }
    fn sub(&self, other: &Self) -> Self {
        unimplemented!()
    }
    fn neg(&self) -> Self {
        Term::new(self.coeff.neg(), self.deg.clone())
    }
    // Ring operations
    fn mul(&self, other: &Self) -> Self {
        Term::new(self.coeff.mul(&other.coeff), self.deg.add(&other.deg))
    }
}

impl<E, P> EuclideanDomain for Term<P> 
    where E: ScalarRing + EuclideanDomain, P: PolyRing<Coeff=E> {

    fn gcd(&self, other: &Self) -> Self {
        Term::new(self.coeff.gcd(&other.coeff), self.deg.gcd(&other.deg))
    }
    fn lcm(&self, other: &Self) -> Self {
        Term::new(self.coeff.lcm(&other.coeff), self.deg.lcm(&other.deg))
    }

    fn divides(&self, other: &Self) -> Option<bool> {
        match (self.coeff.divides(&other.coeff), self.deg.divides(&other.deg)) {
            (Some(a), Some(b)) => Some(a && b),
            _                  => None,
        }
    }
}

impl<F: Field, P: PolyRing<Coeff=F>> Term<P> {
    pub fn div(&self, other: Self) -> Option<Self> {
        if self.divides(&other) == Some(true) {
            Some(Term::new(
                self.coeff.div(&other.coeff).unwrap(),
                self.deg.sub(&other.deg).unwrap()
            ))
        } else {
            None
        }
    }
}

// <><><><><><><><><><> General Polynomial Functions <><><><><><><><><><> //
impl<'a, P: PolyRing> Poly<'a, P> {

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
    pub fn deg(&self)       -> usize     { <P::Var>::tot_deg(&self.lm()) }
    pub fn lm(&self)        -> P::Var    { self.lt().deg }
}

pub trait MonomialOrdering<I: Variate>: Clone + Eq {
    fn cmp(a: &I, b: &I) -> Ordering;
}

pub trait Variate: Zero + Clone + Eq + Debug {
    type NumVar: VarNumber;

    fn tot_deg(&self) -> usize;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Option<Self>;

    // This is implementing the EuclideanDomain trait, but this doesn't actually form a
    // ring so it isn't technically one
    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    fn divides(&self, other: &Self) -> Option<bool>;
}

// Basic Arithmetic operations for polynomial
impl<'a, P: PolyRing> Poly<'a, P> {

    pub fn term_scale(&self, term: &Term<P>) -> Poly<'a, P> {
        let new_terms = self.terms.iter().map(|Term { coeff, deg }| 
            Term::new(
                coeff.mul(&term.coeff),
                deg.add(&term.deg)
            )).collect();
        Poly::from_terms_unchecked(new_terms, self.ring)
        
    }
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
        let (smol, bigg) = match <P::Ord>::cmp(&polya.lm(), &polyb.lm()) {
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
            res.push( match <P::Ord>::cmp(&a.deg, &b.deg) { // Compare their degrees
                Ordering::Less    => { i += 1; a.scale(&smol.ls) },
                Ordering::Greater => { j += 1; b.scale(&bigg.ls) },
                Ordering::Equal   => {
                    i += 1; j += 1;
                    let c = a.coeff.mul(&smol.ls).add(&b.coeff.mul(&bigg.ls));
                    if c == <P::Coeff>::zero() {
                        continue
                    }
                    Term { coeff: c, deg: a.deg.clone() }
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