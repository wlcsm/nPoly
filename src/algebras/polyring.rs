use crate::algebras::*;
use generic_array::*;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::marker::PhantomData;

// <><><><><><><><><><> Poly Ring Domain <><><><><><><><><><> //
pub trait PolyRing: Eq + PartialEq + Clone + std::fmt::Debug {
    type Coeff: ScalarRing;
    type Var: Variate;
    type Ord: MonomialOrdering<Self::Var>;
    fn symb(&self) -> Vec<char>;
}
pub trait VarNumber: ArrayLength<usize> + Eq + PartialEq + Clone + Debug {}

// This is where I use the associated_type_bounds feature
pub trait FPolyRing: PolyRing<Coeff: Field> {}
pub trait PolyRingDiscrete: PolyRing<Coeff: Eq + PartialEq + Ord> {}

use generic_array::typenum::{U2, U3};

impl VarNumber for U2 {}
impl VarNumber for U3 {}

// The reason I have kept this generic for now (when at the moment there isn't any immediate
// need to) is because I plan on implementing quotient rings for which I will need access to generic parameters
// I have used a Vector even though it would probably be better to use a Generic array with the
// exact size that I need because you cant just use a generic size type, it also needs to include
// the type of the data that it is storing. So normally I specify that it is storing usize (the
// type of the indices) but here it needs to store a char. Until generic arrays remove this
// requirement I will use a Vector and contracts to ensure the number of symbols is equal to the
// number of indices.
use std::fmt;
#[derive(Clone)]
pub struct PRDomain<R: ScalarRing, U: Variate, M: MonomialOrdering<U>> {
    pub(crate) vars: Vec<char>,
    ring_parameters: PhantomData<(R, U, M)>,
}

impl<F, U, M> PartialEq for PRDomain<F, U, M> 
where
    F: ScalarRing,
    U: Variate,
    M: MonomialOrdering<U>,
{
    fn eq(&self, other: &Self) -> bool {
        self.vars == other.vars
    }
}

impl<F, U, M> Eq for PRDomain<F, U, M> 
where
    F: ScalarRing,
    U: Variate,
    M: MonomialOrdering<U>
{}

impl<F, U, M> FPolyRing for PRDomain<F, U, M>
where
    F: Field,
    U: Variate,
    M: MonomialOrdering<U>,
{
}

impl<R, U, M> Debug for PRDomain<R, U, M>
where
    R: ScalarRing,
    U: Variate,
    M: MonomialOrdering<U>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vars)
    }
}

impl<R, U, M> PolyRing for PRDomain<R, U, M>
where
    R: ScalarRing,
    U: Variate,
    M: MonomialOrdering<U>,
{
    type Coeff = R;
    type Var = U;
    type Ord = M;
    fn symb(&self) -> Vec<char> {
        self.vars.clone()
    }
}

use generic_array::typenum::Unsigned;

impl<R, U, M> PRDomain<R, U, M>
where
    R: ScalarRing,
    U: Variate,
    M: MonomialOrdering<U>,
{
    pub fn new(vars: Vec<char>) -> Self {
        if vars.len() != U::NumVar::to_usize() {
            panic!("Number of variable symbols supplied to create ring is not the same as its type definition");
        }
        PRDomain {
            vars,
            ring_parameters: PhantomData,
        }
    }
}

// <><><><><><><><><><> Polynomial <><><><><><><><><><> //
#[derive(Debug, Clone)]
pub struct Poly<'a, P: PolyRing> {
    pub(crate) terms: Vec<Term<P>>,
    pub(crate) ring: &'a P,
}

impl<'a, P: PolyRing> Poly<'a, P> {
    // TODO Probably should make a checked version where I sort the list
    // and remove duplicates
    pub fn from_terms(terms: Vec<Term<P>>, ring: &'a P) -> Poly<'a, P> {
        // Sort
        let mut terms_sorted = terms.clone();
        terms_sorted.sort_by(|a, b| <P::Ord>::cmp(&a.deg, &b.deg));
        // Remove duplicates
        let mut term_no_duplicates: Vec<Term<P>> = Vec::new();
        for el in terms_sorted {
            if let Some(last) = term_no_duplicates.last_mut() {
                if el.deg == last.deg {
                    last.coeff = last.coeff + el.coeff;
                    continue;
                }
            }
            term_no_duplicates.push(el)
        }
        // Filter out all zero terms
        let final_terms = term_no_duplicates
            .into_iter()
            .filter(|a| a.coeff != P::Coeff::zero())
            .collect();

        Poly::from_terms_unchecked(final_terms, ring)
    }

    pub fn from_terms_unchecked(terms: Vec<Term<P>>, ring: &'a P) -> Poly<'a, P> {
        Poly {
            terms,
            ring,
        }
    }
}
#[cfg(test)]
mod tests {

    use super::*;
    // use crate::algebras::real::RR;
    // use crate::parse::*;
    // use crate::polyu::*;

    // #[test]
    // fn from_terms_test() {
    //     let ring = PRDomain::<RR, UniIndex, UnivarOrder>::new(vec!['x']);
    //     let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98 - 6.0x^2").unwrap();
    //     println!("{:?}", a);
    //     println!("{}", a);
    // }
}
// <><><><><><><><><><> Term <><><><><><><><><><> //
#[derive(Debug, Clone)]
pub struct Term<P: PolyRing> {
    pub coeff: P::Coeff,
    pub deg: P::Var,
}

use std::cmp::{Ord, PartialOrd};

impl<P: PolyRingDiscrete> PartialEq for Term<P> {
    fn eq(&self, other: &Self) -> bool {
        self.deg == other.deg && self.coeff == other.coeff
    }
}

impl<P: PolyRingDiscrete> Eq for Term<P> {}

impl<P: PolyRingDiscrete> Ord for Term<P> {
    fn cmp(&self, other: &Self) -> Ordering {
        <P::Ord>::cmp(&self.deg, &other.deg)
    }
}

use alga::general::{Identity, Additive, AbstractGroup, Multiplicative};

impl<P: PolyRingDiscrete> PartialOrd for Term<P> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: PolyRing> Term<P> {
    pub fn new(coeff: P::Coeff, deg: P::Var) -> Self {
        Term { coeff, deg }
    }

    // pub fn scale(&self, scalar: &P::Coeff) -> Self {
    //     Term::new(self.coeff.mul(scalar), self.deg.clone())
    // }
    pub fn cmpdeg(&self, other: Self) -> Ordering {
        <P::Ord>::cmp(&self.deg, &other.deg)
    }
}
impl<P: PolyRing> Zero for Term<P> {
    fn zero() -> Self {
        Term::new(<P::Coeff>::zero(), <P::Var>::zero())
    }
}

// impl<P: PolyRing> One for Term<P> {
//     fn one() -> Self {
//         Term::new(<P::Coeff>::one(), <P::Var>::zero())
//     }
// }

impl<T: AbstractGroup<Multiplicative>> One for T {
    fn one() -> Self {
        <Self as Identity<Multiplicative>>::identity()
    }
}

impl<T: AbstractGroup<Additive>> Zero for T {
    fn zero() -> Self {
        <Self as Identity<Additive>>::identity()
    }
}

impl<P: PolyRing> Term<P> {
    fn neg(&self) -> Self {
        Term::new(-self.coeff, self.deg.clone())
    }
    // Ring operations
    fn mul(&self, other: &Self) -> Self {
        Term::new(self.coeff * other.coeff, self.deg.add(&other.deg))
    }
}

impl<E, P> Term<P>
where
    E: ScalarRing + EuclideanDomain,
    P: PolyRing<Coeff = E>,
{
    fn gcd(&self, other: &Self) -> Self {
        Term::new(self.coeff.gcd(&other.coeff), self.deg.gcd(&other.deg))
    }
    fn lcm(&self, other: &Self) -> Self {
        Term::new(self.coeff.lcm(&other.coeff), self.deg.lcm(&other.deg))
    }

    fn divides(&self, other: &Self) -> Option<bool> {
        // Checks if a term is divisible by another term. For this to occur both the monomial and the coefficients, 
        // need to divide each other. This returns a zero if the coefficient of the divisor is zero.
        match (
            self.coeff.divides(&other.coeff),
            self.deg.divides(&other.deg),
        ) {
            (Some(a), Some(b)) => Some(a && b),
            // Note that a index dividing an index will always return a Some() value
            (None, _) => None,
            _ => Some(false),
        }
    }
}

impl<P: PolyRing> Term<P> {
    pub fn to_str(&self, ring: &P) -> String {
        let mut term = self.coeff.to_string();

        for (i, symb) in ring.symb().iter().enumerate() {
            match self.deg.get(i).unwrap() {
                0 => {}
                1 => term.push(*symb),
                d => term.push_str(&format!("{}^{}", symb, d)),
            }
        }
        term
    }
}

impl<P: FPolyRing> Term<P> {
    // Evaluates self / other
    pub fn div(&self, other: Self) -> Option<Self> {
        if other.divides(&self) == Some(true) {
            Some(Term::new(
                self.coeff.div(&other.coeff).unwrap(),
                self.deg.sub(&other.deg).unwrap(),
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
        unsafe { self.terms.get_unchecked(i) }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<P> {
        match self.num_terms() {
            0 => Term::zero(),
            n => {
                self.get(n - 1).unwrap().clone()
            }
        }
    }

    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }
    pub fn lc(&self) -> P::Coeff {
        self.lt().coeff
    }
    pub fn deg(&self) -> usize {
        <P::Var>::tot_deg(&self.lm())
    }
    pub fn lm(&self) -> P::Var {
        self.lt().deg
    }

    pub fn has(&self, t: &P::Var) -> Option<P::Coeff> {
        // Does a binary search for the term, returns the coefficient if
        // it was found and nonzero
        match self.terms.binary_search_by(|a| <P::Ord>::cmp(&a.deg, &t)) {
            Ok(i) => {
                if self.terms[i].coeff != <P::Coeff>::zero() {
                    Some(self.terms[i].coeff)
                } else {
                    None
                }
            }
            Err(_) => None,
        }
    }
}

pub trait MonomialOrdering<I: Variate>: Clone + Eq {
    fn cmp(a: &I, b: &I) -> Ordering;
}

pub trait Variate: Zero + Clone + Eq + Debug + Ord {
    type NumVar: VarNumber;

    fn get(&self, ind: usize) -> Option<&usize>;
    fn set(&mut self, ind: usize, val: usize) -> Option<()>;
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
        let new_terms = self
            .terms
            .iter()
            .map(|Term { coeff, deg }| Term::new(*coeff + term.coeff, deg.add(&term.deg)))
            .collect();
        Poly::from_terms_unchecked(new_terms, self.ring)
    }
    pub fn is_zero(&self) -> bool {
        self.num_terms() == 0
            || (self.num_terms() == 1 && self.get(0).unwrap().coeff == <P::Coeff>::zero())
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
        self.elementwise_map(|t| t.neg())
    }

    pub fn zero(&self) -> Self {
        Poly::from_terms_unchecked(vec![], self.ring)
    }

    pub fn scale(&self, scalar: P::Coeff) -> Self {
        if scalar == <P::Coeff>::zero() {
            Poly::zero(&self)
        } else {
            self.elementwise_map(|t| Term::new(t.coeff * scalar, t.deg.clone()))
        }
    }

    pub fn scale_ass(&mut self, scalar: P::Coeff) {
        if scalar == <P::Coeff>::zero() {
            self.clone_from(&Poly::zero(&self));
        } else {
            self.elementwise_map_mut(|t| t.coeff.mul_ass(&scalar))
        }
    }

    pub fn elementwise_map(&self, map: impl Fn(&Term<P>) -> Term<P>) -> Self {
        Poly::from_terms(self.terms.iter().map(map).collect(), self.ring)
    }

    pub fn elementwise_map_mut(&mut self, map: impl Fn(&mut Term<P>)) {
        self.terms.iter_mut().for_each(map)
    }

    pub fn elementwise_add(polya: &Self, polyb: &Self) -> Self {
        if polya.ring != polyb.ring {
            panic!("Try to add two poly with different rings I don't know how to do that yet")
        }

        // Knowing which polynomial is bigger is advantages here
        let (smol, bigg) = match polya.num_terms().cmp(&polyb.num_terms()) {
            Ordering::Less => (&polya, &polyb),
            _ => (&polyb, &polya),
        };

        // Note: with_capacity allocates memory, should be deallocate some after we
        //       finished?
        let mut res = Vec::with_capacity(bigg.num_terms() + smol.num_terms());

        let mut i = 0;
        let mut j = 0;

        while i < smol.num_terms() {
            let (a, b) = (smol.get_uc(i), bigg.get_uc(j));
            res.push(
                match a.deg.cmp(&b.deg) {
                    // Compare their degrees
                    Ordering::Less => {
                        i += 1;
                        a.clone()
                    }
                    Ordering::Greater => {
                        j += 1;
                        b.clone()
                    }
                    Ordering::Equal => {
                        i += 1;
                        j += 1;
                        let c = a.coeff + b.coeff;
                        if c == <P::Coeff>::zero() {
                            continue;
                        }
                        Term {
                            coeff: c,
                            deg: a.deg.clone(),
                        }
                    }
            })
        }

        // Append any remaining terms to the result vector
        for k in j..bigg.num_terms() {
            res.push(bigg.get_uc(k).clone())
        }

        Poly::from_terms_unchecked(res, polya.ring)
    }
}
