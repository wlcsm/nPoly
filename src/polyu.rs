use crate::algebras::*;
use crate::algebras::polyring::*;
use std::cmp::Ordering;

#[derive(Eq, PartialEq, Clone, Debug)]
pub(crate) struct Univariate(String);

impl<T: ScalarRing> PRDomain<T, Univariate> {
    pub fn univar(symb: String) -> PRDomain<T, Univariate> {
        PRDomain::new(Univariate(symb))
    }
}

impl Variate for Univariate {
    fn cmp(&self, a: &TermIndex, b: &TermIndex) -> Ordering {
        a.0.cmp(&b.0)
    }
    fn tdeg(&self, index: &TermIndex) -> usize {
        index.0 as usize
    }
    fn zero() -> TermIndex { TermIndex(0) }
}

// <><><><><><><><> Constructors <><><><><><><><> //
impl<'a, T: ScalarRing> Poly<'a, T, Univariate> {

    pub fn from_coeff(ring: &'a PRDomain<T, Univariate>, coeffs: Vec<T>) -> Poly<'a, T, Univariate> {
        // Automatically compress the terms argument
        let terms = coeffs.into_iter().enumerate()
                          .filter(|(_, c)| *c != <T>::zero())
                          .map(|(i, c)| Term::new(c, TermIndex(i as u64)))
                          .collect();

        Poly::from_terms_unchecked(terms, ring)
    }
}

impl<'a, T: ScalarRing> PolyMul for Poly<'a, T, Univariate> {
    /// Standard O(n^2) multiplication
    fn mul(&self, other: &Self) -> Self {

        let mut acc = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for Term {coeff: c_a, deg: d_a } in self.terms.iter() {
            for Term {coeff: c_b, deg: d_b } in other.terms.iter() {
                acc[d_a.add(&d_b).0 as usize].add_ass(&c_a.mul(&c_b));
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