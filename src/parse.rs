// Parses a string into a polynomial
// Going to use a trait object to infer the type, though it would be good if we
// could specify the types later in our function.

extern crate regex;

use crate::algebras::polyring::*;
use crate::algebras::*;
use crate::error::PolyErr;
use num_traits::identities::Zero;

use regex::Regex;

// Made my own one because I currently can't use the std::str::FromStr because I need to
// specify what the ring of the polynomial should be beforehand. Since I cannot include the indeterminate's
// symbols at the moment, I can't infer the indeterminates symbols from the type signature
pub trait MyFromStr<'a, P: PolyRing>: Sized {
    type Err;
    fn from_str(ring: &'a P, s: &str) -> Result<Self, Self::Err>;
}



// TODO implement this
// impl<'a, P: PolyRing> MyFromStr<'a, P> for Term<P> {
//     type Err = PolyErr;

//     fn from_str(ring: &'a P, s: &str) -> Result<Self, Self::Err> {

//         let mono = String::new();
//         let mono_named = String::new();

//         // Build the named and unnamed regexes
//         for s in ring.symb() {
//             mono.push_str(&format!(r"({var}\^\d*)?", var = s));
//             mono_named.push_str(&format!(r"(?:{var}\^(?P<{var}>\d*))?", var = s));
//         }

//         if !mono.is_match(s) {
//             Err(PolyErr::ParsePolyError)
//         } else {
//             let mut acc = Vec::new();
//             for caps in mono_named.captures_iter(s) {

//                 // Extract the indices of the indeterminates in the term
//                 let mut degrees = P::Mon::zero();
//                 for (i, symb) in ring.symb().iter().enumerate() {
//                     if let Some(n) = caps.name(symb.to_string().as_str()) {
//                         degrees.set(
//                             i,
//                             n.as_str()
//                                 .parse::<usize>()
//                                 .map_err(|_| PolyErr::ParsePolyError)?,
//                         );
//                     }
//                 }
//                 acc.push(Monomial::new(degrees))
//             }
//             // println!("acc = {:?}", acc);
//             Ok(Poly::<'a, P>::from_terms(acc, Some(ring)))
//         }
//     }
// }
impl<'a, P: PolyRing> MyFromStr<'a, P> for Poly<'a, P> {
    type Err = PolyErr;

    fn from_str(ring: &'a P, s: &str) -> Result<Self, Self::Err> {

        // Build the named and unnamed regexes
        let mut mono = format!("{}", P::Coeff::REGEX);
        let mut mono_named = format!("(?P<coeff>{})", P::Coeff::REGEX);

        for s in ring.symb() {
            mono.push_str(&format!(r"({var}\^\d*)?", var = s));
            mono_named.push_str(&format!(r"(?:{var}\^(?P<{var}>\d*))?", var = s));
        }

        // Validates the entire thing
        let whole_regex = Regex::new(&format!(r"^{m}(\s*(\+|-)\s*{m})*$", m = mono)).unwrap();
        // Parse each term. The "?" is in there for the first term which might not have a + or -,
        // all the other terms are guarenteed to have one because it is checked in the whole_regex
        let term_regex = Regex::new(&format!(r"(?P<sign>\+|-)?\s*{}", mono_named)).unwrap();

        if !whole_regex.is_match(s) {
            Err(PolyErr::ParsePolyError)
        } else {
            let mut acc = Vec::new();
            for caps in term_regex.captures_iter(s) {
                // Extract the coefficient, at the moment, it needs a 1 there
                // println!("caps = {:?}", caps);
                let mut coeff = caps["coeff"]
                    .parse::<P::Coeff>()
                    .map_err(|_| PolyErr::ParsePolyError)?;

                // Take into account if it is a subtraction
                if let Some("-") = caps.name("sign").map(|c| c.as_str()) {
                    coeff = -coeff;
                }

                // Extract the indices of the interminates in the term
                let mut degrees = P::Mon::zero();
                for (i, symb) in ring.symb().iter().enumerate() {
                    if let Some(n) = caps.name(symb.to_string().as_str()) {
                        degrees.set(
                            i,
                            n.as_str()
                                .parse::<usize>()
                                .map_err(|_| PolyErr::ParsePolyError)?,
                        );
                    }
                }
                acc.push(Term::new(coeff, degrees))
            }
            // println!("acc = {:?}", acc);
            Ok(Poly::<'a, P>::from_terms(acc, Some(ring)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::integers::ZZ;
    use crate::algebras::real::RR;
    use crate::polym::*;
    // use crate::polyu::*;
    use generic_array::typenum::U2;

    #[test]
    fn parsing_test() {
        // Univariate
        // let ring = PRDomain::<ZZ, UniIndex, UnivarOrder>::new(vec!['x']);
        // let a = Poly::from_str(&ring, "3x^2 + 5x^98 + 0x^2 + 1x^2 - 1x^2").unwrap();
        // println!("{:?}", a);
        // println!("{}", a);

        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);

        let a = Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap();

        // Multivariate
        let ring = PRDomain::<ZZ, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "3x^2y^6 + 5x^98y^2").unwrap();
        let b = Poly::from_str(&ring, "5x^6").unwrap();
        println!("a = {:?}", a);
        println!("b = {:?}", b);
        println!("a = {}", a);
        println!("b = {}", b);
    }
}
