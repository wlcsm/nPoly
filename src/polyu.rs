use crate::algebras::*;
use crate::algebras::polyring::*;
use std::cmp::Ordering;
use generic_array::arr;

impl<R: ScalarRing> PRDomain<R, UniIndex, UnivarOrder> {
    pub fn univar(symb: usize) -> PRDomain<R, UniIndex, UnivarOrder> {
        PRDomain::new(arr![usize; symb])
    }
}
// Shorthand
pub type PolyU<'a, R> = Poly<'a, PRDomain<R, UniIndex, UnivarOrder>>;

pub trait PolyRingUni: PolyRing<Var=UniIndex, Ord=UnivarOrder> {}

impl<R: ScalarRing> PolyRingUni for PRDomain<R, UniIndex, UnivarOrder> {}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct UnivarOrder();
impl MonomialOrdering<UniIndex> for UnivarOrder {
    fn cmp(a: &UniIndex, b: &UniIndex) -> Ordering {
        a.0.cmp(&b.0)
    }
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub struct UniIndex ( pub(crate) usize );

impl Zero for UniIndex {
    fn zero() -> Self { UniIndex(0) }
}
use std::cmp::{min, max};
use typenum::U1;
impl VarNumber for U1 {}

impl Variate for UniIndex {
    type NumVar = U1;

    fn tot_deg(self: &Self) -> usize {
        self.0
    }

    fn add(&self, other: &Self) -> Self {
        UniIndex(self.0 + other.0)
    }

    fn sub(&self, other: &Self) -> Option<Self> {
        if other.0 <= self.0 {
            Some(UniIndex(self.0 - other.0))
        } else {None}
    }
    
    fn divides(&self, other: &Self) -> Option<bool> {
        Some(self.0 <= other.0)
    }
    fn gcd(&self, other: &Self) -> Self {
        UniIndex(min(self.0, other.0))
    }
    fn lcm(&self, other: &Self) -> Self {
        UniIndex(max(self.0, other.0))
    }
}

// <><><><><><><><> Constructors <><><><><><><><> //
impl<'a, P: PolyRingUni> Poly<'a, P> {

    pub(crate) fn from_coeff(ring: &'a P, coeffs: Vec<P::Coeff>) -> Poly<'a, P> {
        // Automatically compress the terms argument
        let terms = coeffs.into_iter().enumerate()
                          .filter(|(_, c)| *c != <P::Coeff>::zero())
                          .map(|(i, c)| Term::new(c, UniIndex(i)))
                          .collect();

        Poly::from_terms_unchecked(terms, ring)
    }
}

impl<'a, F: Field> PolyU<'a, F> {
    /// Standard O(n^2) multiplication
    pub fn mul(&self, other: &Self) -> Self {

        let mut acc = vec![<F>::zero(); self.deg() + other.deg() + 1];

        for Term {coeff: c_a, deg: d_a } in self.terms.iter() {
            for Term {coeff: c_b, deg: d_b } in other.terms.iter() {
                // I still don't know why I can't use method syntax for the tdeg calls...
                acc[<UniIndex>::tot_deg(d_a) + <UniIndex>::tot_deg(d_b)]
                    .add_ass(&c_a.mul(&c_b));
            }
        }

        Poly::from_coeff(self.ring, acc)
    }
}


// impl fmt::Display for Poly<ZZ> {

//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         let sign = |x: ZZ| { if x < 0 {" - "} else {" + "} };
//         // Formats a nomial: Assumes that num is not zero
//         let nomial = |num: &Monomial::<ZZ>| -> String {format!("{}{}",
//             if num.coeff.abs() == 1 {"".to_string()}  else {format!("{}", num.coeff.abs())},
//             if num.deg         == 1 {"x".to_string()} else {format!("x^{}", num.deg)}
//         )};

//         // Perform an extra check on the first element.
//         // This is only one where degree can be zero.
//         let mut acc: String =
//             if self.terms[0].deg == 0 {
//                 format!("{}", self.terms[0].coeff)
//             } else {
//                 nomial(&self.terms[0])
//             };

//         self.terms.iter().skip(1)
//                          .for_each(|x|
//                                    acc.push_str(&format!("{}{}", sign(x.coeff),
//                                                          nomial(x))));

//         write!(f, "{}", acc)
//     }

// }

// impl std::str::FromStr for Poly<ZZ> {
//     /// The function to parse a string into a polynomial type
//     type Err = PolyErr;

//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         // Clean and remove square brackets
//         let poly_iter = s[1..].trim()
//                          .trim_matches(|p| p == '[' || p == ']' )
//                          .split(',');

//         // Parse each element into i32.
//         let mut acc: Vec<Monomial<ZZ>> = Vec::new();
//         for (i, x) in poly_iter.enumerate() {
//             acc.push(Monomial::<ZZ>::new(ZZ(x.parse::<i32>()?), i))
//         };

//         Ok(Poly::from_monomials(None, acc)?)
//     }
// }