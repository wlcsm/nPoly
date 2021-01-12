use crate::algebras::*;
use generic_array::*;
use std::cmp::Ordering;
use std::fmt;
use std::fmt::Debug;
use std::marker::PhantomData;

// Polynomial Ring Domain
pub trait PolyRing: Eq + PartialEq + Clone + std::fmt::Debug {
    type Coeff: ScalarRing;
    //type Ord: MonOrd<Index = Self::Mon>;
    type Mon: Monomial<NumVar = Self::NumVar>;
    type NumVar: VarNumber;
    fn symb(&self) -> Vec<char>;
}

pub trait VarNumber: ArrayLength<usize> + Eq + PartialEq + Clone + Debug + std::hash::Hash {}

// Super trait for a polynomial whose coefficients are in a field
// This is where I use the associated_type_bounds feature
pub trait FPolyRing: PolyRing<Coeff: ScalarField> {}

// Super trait for a polynomial whose coefficients are discrete
pub trait PolyRingDiscrete: PolyRing<Coeff: Eq + PartialEq + Ord> {}

use generic_array::typenum::{U2, U3, U10};

impl VarNumber for U2 {}
impl VarNumber for U3 {}
impl VarNumber for U10 {}

// I used a vector here instead of an array with exact length because GenericArray
// was a real pain to deal with at the time.
#[derive(Clone)]
pub struct PRDomain<R: ScalarRing, M: MonOrd> {
    pub(crate) vars: Vec<char>,
    ring_parameters: PhantomData<(R, M)>,
}

// Can't use the derive macro it seems because it doesn't know how to compare PhantomData perhaps
impl<F, M> PartialEq for PRDomain<F, M>
where
    F: ScalarRing,
    M: MonOrd,
{
    fn eq(&self, other: &Self) -> bool {
        self.vars == other.vars
    }
}

impl<F, M> Eq for PRDomain<F, M>
where
    F: ScalarRing,
    M: MonOrd,
{
}

impl<F, M> FPolyRing for PRDomain<F, M>
where
    F: ScalarField,
    M: MonOrd,
{
}

impl<R, M> Debug for PRDomain<R, M>
where
    R: ScalarRing,
    M: MonOrd,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vars)
    }
}

impl<R, M> PolyRing for PRDomain<R, M>
where
    R: ScalarRing,
    M: MonOrd,
{
    type Coeff = R;
    type Mon = M::Index;
    type NumVar = <M::Index as Monomial>::NumVar;
    type Ord = M;
    fn symb(&self) -> Vec<char> {
        self.vars.clone()
    }
}

use generic_array::typenum::Unsigned;

impl<R, M> PRDomain<R, M>
where
    R: ScalarRing,
    M: MonOrd,
{
    /// Creates a PRDomain instance from a list of characters denoting the symbols
    /// of the indeterminates.
    /// Panics if the number of characters supplied is not equal to the number of indeterminates
    /// specified in the MonOrd trait.
    /// This could be checked at compile time if we only took GenericArrays instead of vectors
    /// but they are a pain at the moment.
    pub fn new(vars: Vec<char>) -> Self {
        if vars.len() != <M::Index as Monomial>::NumVar::to_usize() {
            panic!("Number of variable symbols supplied to create ring is not the same as its type definition");
        }
        PRDomain {
            vars,
            ring_parameters: PhantomData,
        }
    }
}

// <><><><><><><><><><> Polynomial <><><><><><><><><><> //

#[derive(Debug, Clone, PartialEq)]
pub struct Poly<'a, P: PolyRing> {
    pub(crate) terms: Vec<Term<P>>,
    pub(crate) ring: Option<&'a P>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyDense<'a, P: PolyRing> {
    pub(crate) terms: Vec<P::Coeff>,
    pub(crate) ring: Option<&'a P>,
}


impl<'a, P: PolyRing> Poly<'a, P> {
    /// Returns a polynomial from a vector of terms.
    /// Sorts the terms with respect to a monomial order and removes any zero terms
    /// ~PO We create no_dup and sorted, could look at only operating on the same piece of memory
    /// for space efficiency.
    pub fn from_terms(terms: Vec<Term<P>>, ring: Option<&'a P>) -> Poly<'a, P> {

        // Remove zero terms then sort
        let mut sorted: Vec<Term<P>> = terms.into_iter().filter(|x| !x.is_zero()).collect();
        sorted.sort_by(|a, b| <P::Ord>::cmp(&a.mon, &b.mon));

        // Remove duplicates
        // TODO Test which is faster
//        let mut no_dup: Vec<Term<P>> = Vec::with_capacity(sorted.len());
//        for el in sorted {
//            if let Some(last) = no_dup.last_mut() {
//                if el.mon == last.mon {
//                    last.coeff = last.coeff + el.coeff;
//                    continue;
//                }
//            }
//            no_dup.push(el)
//        }


        let mut last = Term::zero();
        let no_dup = sorted.into_iter().filter_map(|a| 
                                      if a.mon == last.mon { 
                                          last.coeff += a.coeff; 
                                          None 
                                      } else { 
                                          let tmp = last.clone(); 
                                          last = a;
                                          Some(tmp) 
                                      }).collect();

        Poly::from_terms_unchecked(no_dup, ring)

    }

    pub fn from_terms_unchecked(terms: Vec<Term<P>>, ring: Option<&'a P>) -> Poly<'a, P> {
        Poly { terms, ring }
    }
}

/// Evaluate the polynomial at some value
impl<'a, P: PolyRing> PolyDense<'a, P> {
    fn eval(&self, val: P::Coeff) -> P::Coeff {
        self.terms.rev().fold(P::Coeff::zero(), |acc, t| {
            acc * val + t
        })
    }
}

// TODO Not the best
impl<'a, P: PolyRing> MulAssign for Poly<'a, P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs
    }
}

//impl<'a, P: PolyRingUni> PolyU<'a, P> {
//    fn eval(&self, val: P::Coeff) -> P::Coeff {
//
//        let prev = 
//        self.terms.rev().fold(P::Coeff::zero(), |acc, t| {
//            (acc + t.coeff) * val.exp(
//
//
//    }
//}

// <><><><><><><><><><> General Polynomial Functions <><><><><><><><><><> //
impl<'a, P: PolyRing> Poly<'a, P> {
    //
    // A get_unchecked implementation
    pub(crate) fn get_unchecked(&self, i: usize) -> &Term<P> {
        unsafe { self.terms.get_unchecked(i) }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<P> {
        match self.num_terms() {
            0 => Term::zero(),
            n => self.get_unchecked(n - 1).clone(),
        }
    }

    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    pub fn lc(&self) -> P::Coeff {
        self.lt().coeff
    }

    pub fn deg(&self) -> usize {
        <P::Mon>::tot_deg(&self.lm())
    }

    pub fn lm(&self) -> P::Mon {
        self.lt().mon
    }

    /// Does a binary search for the term and returns the coefficient if
    /// it was found and nonzero
    pub fn has(&self, t: &P::Mon) -> Option<P::Coeff> {
        match self.terms.binary_search_by(|a| <P::Ord>::cmp(&a.mon, &t)) {
            Ok(i) => Some(self.terms[i].coeff),
            Err(_) => None,
        }
    }
}



// <><><><><><><><><><> Term <><><><><><><><><><> //
/// The Term struct.
/// It holds a coefficient and a monomial, both are generic.
/// Terms are not ordered but monomials are.
/// This is because two terms should be equal if and only if both their coefficients and
/// monomials are the same, but ordering only depends on the monomials.
///
/// This might not cause problems but I'm hesitant. I think its best to just specify how we want
/// the terms to be ordered during the application, not as an intrinsic property
#[derive(Debug, Clone, PartialEq)]
pub struct Term<P: PolyRing> {
    pub coeff: P::Coeff,
    pub mon: P::Mon,
}

//use std::cmp::Ord;

impl<P: PolyRing> Eq for Term<P> {}

impl<P: PolyRing> Term<P> {
    pub fn new(coeff: P::Coeff, mon: P::Mon) -> Self {
        Term { coeff, mon }
    }
}

use num_traits::{One, Zero};
use std::ops::{Mul, Add, Sub, MulAssign};

impl<P: PolyRing> Mul for Term<P> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Term::new(self.coeff * other.coeff, self.mon + other.mon)
    }
}

impl<P: PolyRing> MulAssign for Term<P> {
    fn mul_assign(&mut self, other: Self) {
        self.coeff *= other.coeff;
        self.mon += other.mon;
    }
}

impl<P: PolyRing> Term<P> {
    pub fn zero() -> Self {
        Term::new(P::Coeff::zero(), P::Mon::zero())
    }

    pub fn is_zero(&self) -> bool {
        self.coeff.is_zero()
    }
}

impl<P: PolyRing> One for Term<P> {
    fn one() -> Self {
        Term::new(P::Coeff::one(), P::Mon::zero())
    }
}


//impl<P: FPolyRing> EuclidDiv for Term<P> {
//    /// Euclidean division "self / other" returning a quotient and remainder.
//    /// Returns None if other is zero
//    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)> {
//        if other.is_zero() {
//            return None
//        } 
//        match self.mon.div(&other.mon) {
//            Some(q) => Some((Term::new(self.coeff / other.coeff, q), Term::zero())),
//            None    => Some((Term::zero(), self.clone())),
//        }
//    }
//}
//
//impl<P: FPolyRing> Term<P> {
//    pub fn gcd(&self, other: &Self) -> Self {
//        Term::new(self.coeff.gcd(&other.coeff), self.mon.gcd(&other.mon))
//    }
//    pub fn lcm(&self, other: &Self) -> Self {
//        Term::new(self.coeff.lcm(&other.coeff), self.mon.lcm(&other.mon))
//    }
//
//    /// Checks if "other divides self" i.e. "other | self".
//    /// ~PO Can be optimised if you really want to
//    pub fn divides(&self, other: &Self) -> Option<bool> {
//        self.euclid_div(other).map(|(_, r)| r.is_zero())
//    }
//}


// The supertraits are here are so that I can use the derive macros
// For some reason even though the MonOrd is going in a PhantomData it
// still doesn't allow me to derive
//pub trait MonOrd: Clone + PartialEq + Debug + Eq + std::hash::Hash {
//    type Index: Monomial;
//    fn cmp(a: &Self::Index, b: &Self::Index) -> Ordering;
//}

pub trait Monomial:
    IntoIterator<Item = usize, IntoIter = GenericArrayIter<usize, <Self as Monomial>::NumVar>>
    + Zero
    + Clone
    + EuclideanDomain
    + Eq
    + Debug
    + std::hash::Hash
    + std::iter::FromIterator<usize>
{
    type NumVar: VarNumber;

    fn get(&self, ind: usize) -> Option<&usize>;
    fn tot_deg(&self) -> usize;
}

/// Implements basic arithmetic operations
#[macro_export]
macro_rules! impl_arithmetic {
    ($op_trait:ident, $op_func:ident, $op_type:ty, $op_symb:expr) => {
        impl<'a, P: PolyRing>  $op_trait for $op_type {
            type Output = Poly<'a, P>;

            fn $op_func(self, other: Self) -> Poly<'a, P> {
                Poly::elementwise_map(&self, &other, |a, b| $op_symb)
            }
        }
    }
}

impl_arithmetic![Add, add, Poly , a + b];
impl_arithmetic![Add, add, &Poly, a + b];
impl_arithmetic![Sub, sub, Poly , a - b];
impl_arithmetic![Sub, sub, &Poly, a - b];

impl<'a, P: PolyRing> Zero for Poly<'a, P> {
    fn zero() -> Self {
        Poly {
            terms: vec![],
            ring: None,
        }
    }
    fn is_zero(&self) -> bool {
        self.num_terms() == 0
    }
}


impl<R: ScalarRing> Mul for PolyDense<R> {
    type Output = Self;

    /// Basic n^2 multiplication algorithm
    fn mul(self, rhs: Self) -> Self::Output {
        let d_a = self.0.len() + 1;
        let d_b = rhs.0.len() + 1;
        let d_res = d_a + d_b - 1;

        let mut res = vec![R::zero(); d_res];

        // We could try to reword this to do all the sums for a particular entry of "res" in one
        // go, rather than iterating over res multiply time. I saw David Harvey do this in Sage
        for i in 0..d_a {
            for j in 0..d_b {
                res[i + j] += self.0[i] * rhs.0[j]
            }
        }

        // Remove any remaining zeros on the end
        while res.last().unwrap().is_zero() {
            res.pop();
        }

        Dense(res)
    }
}


//use std::collections::HashMap;
//
//impl<'a, P: PolyRing> Mul for Poly<'a, P> {
//    type Output = Self;
//
//    fn mul(self, rhs: Self) -> Self::Output {
//        if self.ring != rhs.ring {
//            panic!("Can't multiply polynomials in different rings yet");
//        }
//
//        let mut res_hash = HashMap::<P::Mon, P::Coeff>::new();
//
//        for (a, b) in iproduct!(self.terms, rhs.terms) {
//            // Inserts the value into the hashmap if empty, or adds it to the current value if not
//            let Term { coeff: c, mon: m } = a * b;
//            res_hash.entry(m).and_modify(|v| *v += c).or_insert(c);
//        }
//
//        let mut res_vec: Vec<Term<P>> =
//            res_hash.into_iter().map(|(k, v)| Term::new(v, k)).collect();
//
//        res_vec.sort_by(|a, b| P::Ord::cmp(&a.mon, &b.mon));
//        Poly::from_terms(res_vec, self.ring)
//    }
//}
//


#[cfg(test)]
mod tests {
    extern crate test;
    use super::*;
    use crate::algebras::real::RR;
    use crate::parse::*;
    use crate::polyu::UniVarOrder;
    use crate::bench::*;
    use crate::polym::*;

    #[test]
    fn schoolbook_mult_test() {
        let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "1.0x^1 + 3.0x^2 + 5.0x^5").unwrap();
        let b = Poly::from_str(&ring, "1.0x^1 + 2.0x^3 + 2.0x^5").unwrap();

        let res = a * b;
        println!("{}", res)
    }

    #[test]
    fn schoolbook_mult_bench() {

        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);

        let mut tester = Tester {
            num_el: 5,
            deg_sizes: arr![usize; 2, 3],
            coeff_size: RR(10.0),
            mult_func: Poly::ref_mul,
            ring: &ring,
        };

        let time = Poly::mult_tester(tester.clone());
        println!("Hashmap time = {} microseconds", time.num_microseconds().unwrap());
    }
}

impl<'a, P: PolyRing> One for Poly<'a, P> {
    fn one() -> Self {
        Poly::from_terms(vec![Term::new(P::Coeff::one(), P::Mon::zero())], None)
    }
}

//impl<'a, P: PolyRing> MyAddMonoid for Poly<'a, P> {
//    fn ref_add(&self, other: &Self) -> Self {
//        Poly::elementwise_map(&self, &other, |a, b| a + b)
//    }
//}
//
//impl<'a, P: PolyRing> MyAddGroup for Poly<'a, P> {
//    fn ref_sub(&self, other: &Self) -> Self {
//        Poly::elementwise_map(&self, &other, |a, b| a - b)
//    }
//}

// FIXME Needs optimisation
//impl<'a, P: PolyRing> MyMulMonoid for Poly<'a, P> {
//    fn ref_mul(&self, other: &Self) -> Self {
//        self.clone() * other.clone()
//    }
//}

// Multiplication by a single term
impl<'a, P: PolyRing> Mul<Term<P>> for &Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn mul(self, term: Term<P>) -> Poly<'a, P> {
        let new_terms = self
            .terms
            .iter()
            .map(|t| t.clone() * term.clone())
            .collect();
        Poly::from_terms_unchecked(new_terms, self.ring)
    }
}

/// Scale assign
/// FIXME This can have a problem if P::Coeff is not an integral domain and then some of the
/// coefficients can become zero but we don't clean them up here
//impl<'a, P: PolyRing> MulAssign<P::Coeff> for Poly<'a, P> {
//    fn mul_assign(&mut self, scalar: P::Coeff) {
//        if scalar.is_zero() {
//            self.terms.clear();
//        } else {
//            self.elementwise_map_mut(|t| t.coeff *= scalar)
//        }
//    }
//}

// Returns the ring when performing operations between two polynomials.
// If one has a None type then this indicates that it is a constant value and
// so it is given the ring of the other polynomial if possible.
// ~PO In the add function, we could just chose the ring of the polynomials which has the
// largest number of elements, since a None type indicates it is a constant and so
// it only has one element.
// This doesn't work for the case where the polynomials have two different rings.
fn resolve_for_constants<'a, P: PolyRing>(
    polya: &Poly<'a, P>,
    polyb: &Poly<'a, P>,
) -> Option<&'a P> {
    match (polya.ring, polyb.ring) {
        (Some(r1), Some(r2)) => {
            if r1 != r2 {
                panic!("Can't handle different rings yet")
            } else {
                Some(r1)
            }
        }
        (Some(r1), None) => Some(r1),
        (None, Some(r2)) => Some(r2),
        (None, None) => None,
    }
}

impl<'a, P: PolyRing> Poly<'a, P> {
    /// Maps a function across all non-zero coefficients of the polynomial.
    /// ~PO
    pub fn elementwise_map(&self, func: impl Fn(&Term<P>) -> Term<P>) -> Self {
        Poly::from_terms(self.terms.iter().map(func).collect(), self.ring)
    }

    pub fn elementwise_map_mut(&mut self, func: impl Fn(&mut Term<P>)) {
        self.terms.iter_mut().for_each(func)
    }

    /// Adds two polynomials together, but can specify the function used e.g. can also do
    /// subtraction.
    /// Note: It cannot do multiplication. This is because it assumes that f(a, 0) = a, f(0, b) = b
    /// which is false in general for multiplication
    /// At the moment it assumes commutativity.
    /// ~PO Add a second function argument that tells it what to do if one of the terms is zero,
    /// is way pointwise multiplications and subtractions (i.e. pass the .neg() operation) can be
    /// done faster.
    pub fn elementwise_map(
        polya: &Self,
        polyb: &Self,
        op: fn(P::Coeff, P::Coeff) -> P::Coeff,
    ) -> Self {
        // Gets the ring the result is going to be in
        let ring = resolve_for_constants(&polya, &polyb);

        // Note: with_capacity allocates memory, should be deallocate some after we
        //       finished?
        let mut res = Vec::with_capacity(polya.num_terms() + polyb.num_terms());

        let mut i = 0;
        let mut j = 0;

        while i < polya.num_terms() && j < polyb.num_terms() {
            let (a, b) = (polya.get_unchecked(i), polyb.get_unchecked(j));
            res.push(match <P::Ord>::cmp(&a.mon, &b.mon) {
                // Compare their degrees
                Ordering::Less => {
                    i += 1;
                    Term::new(op(a.coeff, P::Coeff::zero()), a.mon.clone())
                }
                Ordering::Greater => {
                    j += 1;
                    Term::new(op(P::Coeff::zero(), b.coeff), b.mon.clone())
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                    let c = op(a.coeff, b.coeff);
                    if c.is_zero() {
                        continue;
                    }
                    Term::new(c, a.mon.clone())
                }
            })
        }

        // Append any remaining terms to the result vector
        for k in i..polya.num_terms() {
            res.push(polya.get_unchecked(k).clone())
        }
        // Append any remaining terms to the result vector
        for k in j..polyb.num_terms() {
            res.push(polyb.get_unchecked(k).clone())
        }

        Poly::from_terms_unchecked(res, ring)
    }
}
