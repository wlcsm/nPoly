use crate::algebras::ScalarRing;
use crate::polydense::PolyDense;

/// Sparse Polynomial
#[derive(Debug, Clone, PartialEq)]
pub struct Poly<R: ScalarRing, const VARS: usize> {
    pub(crate) terms: Vec<Term<R, VARS>>,
}

/// Returns a polynomial from a vector of terms.
/// Sorts the terms with respect to a monomial order and removes any zero terms
impl<R: ScalarRing, const VARS: usize> Poly<R, VARS> {
    pub fn from_terms(terms: Vec<Term<R, VARS>>) -> Poly<R, VARS> {
        // Remove zero terms then sort
        let mut sorted: Vec<Term<R, VARS>> = terms.into_iter().filter(|x| !x.is_zero()).collect();
        sorted.sort_by(|a, b| a.mon.cmp(&b.mon));

        // Remove duplicates
        let mut no_dup: Vec<Term<R, VARS>> = Vec::with_capacity(sorted.len());
        for el in sorted {
            if let Some(last) = no_dup.last_mut() {
                if el.mon == last.mon {
                    last.coeff += el.coeff;
                    continue;
                }
            }
            no_dup.push(el)
        }

        Poly { terms: no_dup }
    }
}

/// Converts a polynomial in sparse representation into its equivalent dense representation
/// by padding with zeros. 
/// If the "deg" parameter has a value, it will pad the vector of coefficients with zero until
/// that number
pub fn to_dense<R: ScalarRing>(input: &Poly<R, 1>, deg: Option<usize>) -> PolyDense<R> {
    let mut result: Vec<R> = Vec::with_capacity(input.lt().mon[0]);

    for Term { coeff, mon } in input.terms.iter() {
        result.resize(mon[0], R::zero());
        result.push(*coeff);
    }

    if let Some(padding) = deg {
        result.resize(padding, R::zero());
    }

    PolyDense { terms: result }
}

impl<R: ScalarRing> Poly<R, 1> {

    /// Creates a sparse polynomial from the input slice
    pub fn from_coeff(coeffs: &[R]) -> Poly<R, 1> {
        let terms = coeffs.into_iter().enumerate().filter_map(|(i, c)| {
            if c.is_zero() {
                None
            } else {
                Some(Term::new(*c, [i]))
            }
        }).collect();
        Poly { terms }
    }

    /// Evaluation of polynomials at a value using Horner's method
    pub fn eval(&self, val: R) -> R {
        let (res, last) = self.terms.iter().rev().fold(
            (R::zero(), self.terms.last().unwrap().mon),
            |(eval, last_mon), t| {
                (
                    eval * num::pow::pow(val, last_mon[0] - t.mon[0]) + t.coeff,
                    t.mon,
                )
            },
        );

        res * num::pow::pow(val, last[0])
    }
}

/// Multiplies sparse polynomials using a HashMap.
/// Though this algorithm is O(n log n), there are many algorithms that may be better in practice.
/// Though batch multiplying many polynomials is quite efficient using this method
impl<R: ScalarRing, const VARS: usize> Mul for Poly<R, VARS> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res_hash: HashMap<[usize; VARS], R> = HashMap::new();

        for (a, b) in iproduct!(self.terms, rhs.terms) {
            // Inserts the value into the hashmap if empty, or adds it to the current value if not
            let Term { coeff: c, mon: m } = a * b;
            res_hash.entry(m).and_modify(|v| *v += c).or_insert(c);
        }

        let mut res_vec: Vec<Term<R, VARS>> =
            res_hash.into_iter().map(|(k, v)| Term::new(v, k)).collect();

        res_vec.sort_by(|a, b| a.mon.cmp(&b.mon));
        Poly::from_terms(res_vec)
    }
}

/// Scaling a polynomial via a term
/// NOTE: May break if isn't an integral domain
impl<R: ScalarRing, const VARS: usize> Mul<Term<R, VARS>> for &Poly<R, VARS> {
    type Output = Poly<R, VARS>;

    fn mul(self, term: Term<R, VARS>) -> Poly<R, VARS> {
        let new_terms = self
            .terms
            .iter()
            .map(|t| t.clone() * term.clone())
            .collect();
        Poly { terms: new_terms }
    }
}

/// Scales a polynomial via a scalar value
/// NOTE: May break if R isn't an integral domain
impl<R: ScalarRing, const VARS: usize> MulAssign<R> for Poly<R, VARS> {
    fn mul_assign(&mut self, scalar: R) {
        if scalar.is_zero() {
            self.terms.clear();
        } else {
            self.terms.iter_mut().map(|t| t.coeff *= scalar).collect()
        }
    }
}

impl<R: ScalarRing, const VARS: usize> MulAssign for Poly<R, VARS> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs
    }
}


impl<R: ScalarRing, const VARS: usize> Poly<R, VARS> {

    pub(crate) fn get_unchecked(&self, i: usize) -> &Term<R, VARS> {
        unsafe { self.terms.get_unchecked(i) }
    }

    // Assumes the lead term is the last element of the vector
    pub fn lt(&self) -> Term<R, VARS> {
        match self.num_terms() {
            0 => Term::zero(),
            n => self.get_unchecked(n - 1).clone(),
        }
    }

    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    pub fn tot_deg(&self) -> usize {
        self.lt().mon.iter().sum()
    }

    /// Does a binary search for the term and returns the coefficient if it was found and nonzero
    pub fn has(&self, t: &[usize; VARS]) -> Option<R> {
        match self.terms.binary_search_by(|a| a.mon.cmp(&t)) {
            Ok(i) => Some(self.terms[i].coeff),
            Err(_) => None,
        }
    }
}

/// The Term struct, consisting on a coefficient and a monomial
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Term<R: ScalarRing, const VARS: usize> {
    pub coeff: R,
    pub mon: [usize; VARS],
}

impl<R: ScalarRing, const VARS: usize> Term<R, VARS> {
    pub fn new(coeff: R, mon: [usize; VARS]) -> Self {
        Term { coeff, mon }
    }
}

use num_traits::{One, Zero};
use std::ops::{Add, Mul, MulAssign, Sub};

impl<R: ScalarRing, const VARS: usize> Mul for Term<R, VARS> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut mon = [0; VARS];
        for i in 0..VARS {
            mon[i] = self.mon[i] + other.mon[i]
        }
        Term::new(self.coeff * other.coeff, mon)
    }
}

impl<R: ScalarRing, const VARS: usize> MulAssign for Term<R, VARS> {
    fn mul_assign(&mut self, other: Self) {
        self.coeff *= other.coeff;
        for i in 0..VARS {
            self.mon[i] += other.mon[i]
        }
    }
}

impl<R: ScalarRing, const VARS: usize> Term<R, VARS> {
    fn zero() -> Self {
        Term::new(R::zero(), [0; VARS])
    }

    fn is_zero(&self) -> bool {
        self.coeff.is_zero()
    }
}

impl<R: ScalarRing, const VARS: usize> One for Term<R, VARS> {
    fn one() -> Self {
        Term::new(R::one(), [0; VARS])
    }
}

impl<R: ScalarRing, const VARS: usize> Zero for Poly<R, VARS> {
    fn zero() -> Self {
        Poly { terms: vec![] }
    }
    fn is_zero(&self) -> bool {
        self.num_terms().is_zero()
    }
}

impl<R: ScalarRing, const VARS: usize> One for Poly<R, VARS> {
    fn one() -> Self {
        Poly::from_terms(vec![Term::new(R::one(), [0; VARS])])
    }
}

/// Macro for implementing basic arithmetic operations
#[macro_export]
macro_rules! impl_arithmetic {
    ($op_trait:ident, $op_func:ident, $op_type:ty) => {
        impl<R: ScalarRing, const VARS: usize> $op_trait for $op_type {
            type Output = Poly<R, VARS>;

            fn $op_func(self, other: Self) -> Poly<R, VARS> {
                Poly::merge(&self, &other, |a, b| a.$op_func(b))
            }
        }
    };
}

impl_arithmetic![Add, add,  Poly<R, VARS>];
impl_arithmetic![Add, add, &Poly<R, VARS>];
impl_arithmetic![Sub, sub,  Poly<R, VARS>];
impl_arithmetic![Sub, sub, &Poly<R, VARS>];

use std::collections::HashMap;

extern crate itertools;
use itertools::iproduct;

use std::cmp::Ordering;

impl<R: ScalarRing, const VARS: usize> Poly<R, VARS> {
    /// Maps a function across all non-zero coefficients of the polynomial.
    pub fn elementwise_map(&self, func: impl Fn(&Term<R, VARS>) -> Term<R, VARS>) -> Self {
        Poly {
            terms: self.terms.iter().map(func).collect(),
        }
    }

    /// Adds two polynomials together, but can specify the function used e.g. can also do
    /// subtraction.
    pub fn merge(polya: &Self, polyb: &Self, op: fn(R, R) -> R) -> Self {
        // Note: with_capacity allocates memory, should be deallocate some after we finished?
        let mut res = Vec::with_capacity(polya.num_terms() + polyb.num_terms());

        let mut i = 0;
        let mut j = 0;

        while i < polya.num_terms() && j < polyb.num_terms() {
            let (a, b) = (polya.get_unchecked(i), polyb.get_unchecked(j));
            let (coeff, mon) = match a.mon.cmp(&b.mon) {
                Ordering::Less => {
                    i += 1;
                    (op(a.coeff, R::zero()), a.mon.clone())
                }
                Ordering::Greater => {
                    j += 1;
                    (op(R::zero(), b.coeff), b.mon.clone())
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                    let c = op(a.coeff, b.coeff);
                    if c.is_zero() {
                        continue;
                    }
                    (c, a.mon.clone())
                }
            };
            res.push(Term::new(coeff, mon))
        }

        // Append any remaining terms to the result vector
        for k in i..polya.num_terms() {
            res.push(polya.get_unchecked(k).clone())
        }
        for k in j..polyb.num_terms() {
            res.push(polyb.get_unchecked(k).clone())
        }

        Poly { terms: res }
    }
}
