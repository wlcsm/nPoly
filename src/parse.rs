// Parses a string into a polynomial
// Going to use a trait object to infer the type, though it would be good if we
// could specify the types later in our function.

extern crate regex;

use crate::algebras::polyring::*;
use crate::algebras::*;
use crate::error::PolyErr;

use regex::Regex;

// Made my own one because I currently can't use the std::str::FromStr because I need to
// specify what the ring of the polynomial should be beforehand. Since I cannot include the indeterminate's
// symbols at the moment, I can't infer the indeterminates symbols from the type signature
pub trait MyFromStr<P: PolyRing>: Sized {
    type Err;
    fn from_str(ring: P, s: &str) -> Result<Self, Self::Err>;
}

impl<P: PolyRing> MyFromStr<P> for Poly<P> {
    type Err = PolyErr;

    fn from_str(ring:  P, s: &str) -> Result<Self, Self::Err> {
        // Build the named and unnamed regexes
        let mut mono = format!("{}", P::Coeff::REGEX);
        let mut mono_named = format!("(?P<coeff>{})", P::Coeff::REGEX);

        for s in ring.symb().unwrap() {
            mono.push_str(&format!(r"({var}\^\d*)?", var = s));
            mono_named.push_str(&format!(r"(?:{var}\^(?P<{var}>\d*))?", var = s));
        }

        // Validates the entire thing
        let whole_regex = Regex::new(&format!(r"^{m}(\s*(\+|-)\s*{m})*$", m = mono)).unwrap();
        // Parse each term. The "?" is in there for the first term which might not have a + or -,
        // all the other terms are guaranteed to have one because it is checked in the whole_regex
        let term_regex = Regex::new(&format!(r"(?P<sign>\+|-)?\s*{}", mono_named)).unwrap();

        if !whole_regex.is_match(s) {
            Err(PolyErr::ParsePolyError)
        } else {
            let mut acc = Vec::new();
            for caps in term_regex.captures_iter(s) {
                // Extract the coefficient, at the moment, it needs a 1 there
                let mut coeff = caps["coeff"]
                    .parse::<P::Coeff>()
                    .map_err(|_| PolyErr::ParsePolyError)?;

                // Take into account if it is a subtraction
                if let Some("-") = caps.name("sign").map(|c| c.as_str()) {
                    coeff = coeff.neg();
                }

                // Extract the indices of the interminates in the term
                let mut degrees = P::Var::zero();
                for (i, symb) in ring.symb().unwrap().iter().enumerate() {
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
            Ok(Poly::<P>::from_terms(acc, ring))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::integers::ZZ;
    use crate::polym::*;
    use crate::polyu::*;
    use generic_array::typenum::U2;

    #[test]
    fn parsing_test() {
        // Univariate
        let ring = PRDomain::<ZZ, UniIndex, UnivarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3x^2 + 5x^98 + 0x^2 + 1x^2 - 1x^2").unwrap();
        println!("{:?}", a);
        println!("{}", a);

        // Multivariate
        let ring = PRDomain::<ZZ, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "3x^2y^6 + 5x^98y^2").unwrap();
        let b = Poly::from_str(&ring, "5x^6").unwrap();
        println!("a = {}", a);
        println!("b = {}", b);
    }
}
