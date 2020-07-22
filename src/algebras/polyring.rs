use crate::algebras::*;
use generic_array::*;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::marker::PhantomData;

// <><><><><><><><><><> Poly Ring Domain <><><><><><><><><><> //
pub trait PolyRing: Eq + PartialEq + Clone + std::fmt::Debug {
    type Coeff: ScalarRing;
    type Mon: Monomial;
    type Ord: MonOrd<Index = Self::Mon>;
    fn symb(&self) -> Vec<char>;
    const ZERO: Self;
}
pub trait VarNumber: ArrayLength<usize> + Eq + PartialEq + Clone + Debug + std::hash::Hash {}

// This is where I use the associated_type_bounds feature
pub trait FPolyRing: PolyRing<Coeff: ScalarField> {}
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
pub struct PRDomain<R: ScalarRing,  M: MonOrd> {
    pub(crate) vars: Vec<char>,
    ring_parameters: PhantomData<(R, M)>,
}

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
    type Ord = M;
    fn symb(&self) -> Vec<char> {
        self.vars.clone()
    }
    /// My solution to the problem of zero and one functions
    const ZERO: PRDomain<R, M> = PRDomain {
        vars: vec![],
        ring_parameters: PhantomData,
    };
}

use generic_array::typenum::Unsigned;

impl<R, M> PRDomain<R, M>
where
    R: ScalarRing,
    M: MonOrd,
{
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

impl<'a, P: PolyRing> Poly<'a, P> {
    /// Returns a polynomial from a vector of terms.
    /// Sorts the terms with respect to a monomial order and removes any zero terms
    /// ~PO
    pub fn from_terms(terms: Vec<Term<P>>, ring: Option<&'a P>) -> Poly<'a, P> {
        // Remove zero terms then sort
        let mut sorted: Vec<Term<P>> = terms.into_iter().filter(|x| !x.is_zero()).collect();
        sorted.sort_by(|a, b| <P::Ord>::cmp(&a.mon, &b.mon));


        // Adds terms with the same degree together
        // acc is an accumulator for all the terms with the same degree, it is then pushed into 
        // the no duplicate vector when it has found all the terms with the same degree
        let mut no_dup: Vec<Term<P>> = Vec::new();
        let mut acc = Term::zero();
        for el in sorted {
            if el.mon == acc.mon {
                acc.coeff += el.coeff;
            } else {
                no_dup.push(acc);
                acc = el;
            }
        }

        // Remove duplicates: This one just looked a little but more messy, I think by not using an 
        // accumulator we might inccur some performance penalties
        // let mut term_no_duplicates: Vec<Term<P>> = Vec::with_capacity(sorted.len());
        // for el in terms_sorted {
        //     if let Some(last) = term_no_duplicates.last_mut() {
        //         if el.mon == last.mon {
        //             last.coeff = last.coeff + el.coeff;
        //             continue;
        //         }
        //     }
        //     term_no_duplicates.push(el)
        // }

        Poly::from_terms_unchecked(no_dup, ring)
    }

    pub fn from_terms_unchecked(terms: Vec<Term<P>>, ring: Option<&'a P>) -> Poly<'a, P> {
        Poly { terms, ring }
    }
}
// #[cfg(test)]
// mod tests {

// use super::*;
// use crate::algebras::real::RR;
// use crate::parse::*;
// use crate::polyu::*;

// #[test]
// fn from_terms_test() {
//     let ring = PRDomain::<RR, MniIndex, UnivarOrder>::new(vec!['x']);
//     let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98 - 6.0x^2").unwrap();
//     println!("{:?}", a);
//     println!("{}", a);
// }
// }
// <><><><><><><><><><> Term <><><><><><><><><><> //
#[derive(Debug, Clone, PartialEq)]
pub struct Term<P: PolyRing> {
    pub coeff: P::Coeff,
    pub mon  : P::Mon,
}

use std::cmp::{Ord, PartialOrd};

// impl<P: PolyRingDiscrete> PartialEq for Term<P> {
//     fn eq(&self, other: &Self) -> bool {
//         self.mon == other.mon && self.coeff == other.coeff
//     }
// }

impl<P: PolyRingDiscrete> Eq for Term<P> {}

impl<P: PolyRingDiscrete> Ord for Term<P> {
    fn cmp(&self, other: &Self) -> Ordering {
        <P::Ord>::cmp(&self.mon, &other.mon)
    }
}

// use alga::general::{Identity, Additive, AbstractGroup, Multiplicative};

impl<P: PolyRingDiscrete> PartialOrd for Term<P> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: PolyRing> Term<P> {
    // FIXME This is somewhat redundant
    pub fn new(coeff: P::Coeff, mon: P::Mon) -> Self {
        Term { coeff, mon }
    }

    // pub fn scale(&self, scalar: &P::Coeff) -> Self {
    //     Term::new(self.coeff.mul(scalar), self.mon.clone())
    // }
    pub fn cmpdeg(&self, other: Self) -> Ordering {
        <P::Ord>::cmp(&self.mon, &other.mon)
    }
}

use num_traits::{One, Zero};

// use alga::general::AbstractMagma;

impl<P: PolyRing> Mul for Term<P> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Term::new(self.coeff * other.coeff, self.mon.ref_add(&other.mon))
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

impl<P: PolyRing> Neg for Term<P> {
    type Output = Self;

    fn neg(self) -> Self {
        Term::new(-self.coeff, self.mon.clone())
    }
}

impl<P: PolyRing> Neg for &Term<P> {
    type Output = Term<P>;

    fn neg(self) -> Term<P> {
        Term::new(-self.coeff, self.mon.clone())
    }
}

impl<P: PolyRing> One for Term<P> {
    fn one() -> Self {
        Term::new(P::Coeff::one(), P::Mon::zero())
    }
}

impl<P: PolyRing> ClosedMul for Term<P> {}

impl<P: PolyRing> MyMulMonoid for Term<P> {
    fn ref_mul(&self, other: &Self) -> Self {
        Term::new(self.coeff * other.coeff, self.mon.ref_add(&other.mon))
    }
}

impl<P: FPolyRing> EuclidDiv for Term<P> {
    /// Euclidean division "self / other" returning a quotient and remainder.
    /// Returns None if other is zero
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)> {
        if other.is_zero() {
            None
        } else {
            match self.mon.div(&other.mon) {
                Some(q) => Some((Term::new(self.coeff / other.coeff, q), Term::zero())),
                None => Some((Term::zero(), self.clone())),
            }
        }
    }
}

impl<P: FPolyRing> Term<P> {
    pub fn gcd(&self, other: &Self) -> Self {
        Term::new(self.coeff.gcd(&other.coeff), self.mon.gcd(&other.mon))
    }
    pub fn lcm(&self, other: &Self) -> Self {
        Term::new(self.coeff.lcm(&other.coeff), self.mon.lcm(&other.mon))
    }

    /// Checks if "self divides other" i.e. "self | other".
    /// Does so by calling the division algorithm can checking if the
    /// remainder is zero, hence it could be optimised
    /// ~PO
    pub fn divides(&self, other: &Self) -> Option<bool> {
        // Checks self divides other.
        self.euclid_div(other).map(|(_, r)| r.is_zero())
    }
}

impl<P: PolyRing> Term<P> {
    pub fn to_str(&self, ring: &Option<&P>) -> String {
        let mut term = self.coeff.to_string();

        // Add the extra variables if it is in a defined ring
        if let Some(r) = ring {
            for (i, symb) in r.symb().iter().enumerate() {
                match self.mon.get(i).unwrap() {
                    0 => {}
                    1 => term.push(*symb),
                    d => term.push_str(&format!("{}^{}", symb, d)),
                }
            }
        }
        term
    }
}

// <><><><><><><><><><> General Polynomial Functions <><><><><><><><><><> //
impl<'a, P: PolyRing> Poly<'a, P> {
    pub(crate) fn get(&self, i: usize) -> Option<&Term<P>> {
        self.terms.get(i)
    }

    // A get_unchecked implementation
    pub(crate) fn get_uc(&self, i: usize) -> &Term<P> {
        unsafe { self.terms.get_unchecked(i) }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<P> {
        match self.num_terms() {
            0 => Term::new(<P::Coeff>::zero(), <P::Mon>::zero()),
            n => self.get(n - 1).unwrap().clone(),
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
                // I feel like we don't need this check anymore since it should be an invariant
                // that there are no zero coefficients 
            // {
                // if self.terms[i].coeff.is_zero() {
                //     Some(self.terms[i].coeff)
                // } else {
                //     None
                // }
            // }
            Err(_) => None,
        }
    }
}

// The supertraits are here are so that I can use the derive macros
// For some reason even though the MonOrd is going in a PhantomData it
// still doesn't allow me to derive
pub trait MonOrd: Clone + PartialEq + Debug {
    type Index: Monomial;
    fn cmp(a: &Self::Index, b: &Self::Index) -> Ordering;
}

// pub trait Indices: MyAddMonoid {
//     type NumVar:  VarNumber;
//     const NUMVAR: usize;
    

//     fn div(&self, other: &Self) -> Option<Self>;

//     fn divides(&self, other: &Self) -> Option<bool> {
//         if self.is_zero() {
//             None
//         } else {
//             Some(self.div(other).is_some())
//         }
//     }
// }

pub trait Monomial: MyAddMonoid + Eq + Debug + std::hash::Hash {
    type NumVar: VarNumber;

    fn get(&self, ind: usize) -> Option<&usize>;
    fn set(&mut self, ind: usize, val: usize) -> Option<()>;
    fn tot_deg(&self) -> usize;

    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;

    fn lex(&self, other: &Self) -> Ordering;

    fn div(&self, other: &Self) -> Option<Self>;

    fn divides(&self, other: &Self) -> Option<bool> {
        if self.is_zero() {
            None
        } else {
            Some(self.div(other).is_some())
        }
    }
    // FIXME This returns None it is isn't divisble and if other is zero
    // there should be something to distinguish these two cases
}

impl<'a, P: PolyRing> Add for &Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn add(self, other: Self) -> Poly<'a, P> {
        Poly::elementwise_add(self, other, |a, b| a + b)
    }
}

impl<'a, P: PolyRing> Add for Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn add(self, other: Self) -> Poly<'a, P> {
        Poly::elementwise_add(&self, &other, |a, b| a + b)
    }
}

impl<'a, P: PolyRing> Sub for &Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn sub(self, other: Self) -> Poly<'a, P> {
        Poly::elementwise_add(self, other, |a, b| a - b)
    }
}

impl<'a, P: PolyRing> Sub for Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn sub(self, other: Self) -> Poly<'a, P> {
        Poly::elementwise_add(&self, &other, |a, b| a - b)
    }
}

impl<'a, P: PolyRing> Neg for &Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn neg(self) -> Poly<'a, P> {
        self.elementwise_map(|t| t.neg())
    }
}

impl<'a, P: PolyRing> Neg for Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn neg(self) -> Poly<'a, P> {
        self.elementwise_map(|t| t.neg())
    }
}

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

impl<'a, P: PolyRing> Mul for Poly<'a, P> {
    type Output = Self;

    fn mul(self, _rhs: Self) -> Self::Output {
        unimplemented!()
    }
}

impl<'a, P: PolyRing> One for Poly<'a, P> {
    fn one() -> Self {
        Poly::from_terms(vec![Term::new(P::Coeff::one(), P::Mon::zero())], None)
    }
}

impl<'a, P: PolyRing> ClosedAdd for Poly<'a, P> {}

impl<'a, P: PolyRing> MyAddMonoid for Poly<'a, P> {
    fn ref_add(&self, other: &Self) -> Self {
        Poly::elementwise_add(&self, &other, |a, b| a + b)
    }
}

impl<'a, P: PolyRing> MyAddGroup for Poly<'a, P> {
    fn ref_sub(&self, other: &Self) -> Self {
        Poly::elementwise_add(&self, &other, |a, b| a - b)
    }
}

impl<'a, P: PolyRing> ClosedMul for Poly<'a, P> {}
// FIXME Needs optimisation
impl<'a, P: PolyRing> MyMulMonoid for Poly<'a, P> {
    fn ref_mul(&self, other: &Self) -> Self {
        self.clone() * other.clone()
    }
}

// Multiplication by a single term
impl<'a, P: PolyRing> Mul<Term<P>> for &Poly<'a, P> {
    type Output = Poly<'a, P>;

    fn mul(self, term: Term<P>) -> Poly<'a, P> {
        let new_terms = self
            .terms
            .iter()
            .map(|Term { coeff, mon }| Term::new(*coeff * term.coeff, mon.ref_add(&term.mon)))
            .collect();
        Poly::from_terms_unchecked(new_terms, self.ring)
    }
}

// Causes some confliction
// impl<'a, P: PolyRing> Mul<P::Coeff> for &Poly<'a, P> {

//     type Output = Poly<'a, P>;

//     fn mul(&self, scalar: P::Coeff) -> Poly<'a, P> {
//         if scalar.is_zero() {
//             Poly::zero()
//         } else {
//             self.elementwise_map(|t| Term::new(t.coeff * scalar, t.mon.clone()))
//         }
//     }
// }

/// Scale assign
/// FIXME This can have a problem if P::Coeff is not an integral domain and then some of the 
/// coefficients can become zero but we don't clean them up here
impl<'a, P: PolyRing> MulAssign<P::Coeff> for Poly<'a, P> {
    fn mul_assign(&mut self, scalar: P::Coeff) {
        if scalar.is_zero() {
            self.terms.clear();
        } else {
            self.elementwise_map_mut(|t| t.coeff *= scalar)
        }
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
    pub fn elementwise_add(polya: &Self, polyb: &Self, func: fn(P::Coeff, P::Coeff) -> P::Coeff) -> Self {
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
            res.push(match <P::Ord>::cmp(&a.mon, &b.mon) {
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
                    let c = func(a.coeff, b.coeff);
                    if c.is_zero() {
                        continue;
                    }
                    Term::new(c, a.mon.clone())
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
