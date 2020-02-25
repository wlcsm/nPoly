use crate::error::PolyErr;
use crate::algebras::*;
use std::marker::PhantomData;



// This bad boi restricts the different possible scalar types
// TODO I can look at implementing so functionality into the coeff trait.
// For example 
// impl<T: PolyRing> Coeff for PolyType {
//     fn add(T, T) -> T;
// }
// I..think that will work?
// Also implementing different Display commands, would be good, for instance we could
// have a function fn fmt_nomial() -> String, which would, in situations such as 
// x(1 + y), carry the x through and return x + xy.

pub enum RingType {Poly, Scalar}

pub trait Coeff {
    fn is_scalar() -> bool;
    fn is_poly() -> bool { !<Self>::is_scalar() }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ScalarType();
impl Coeff for ScalarType {
    fn is_scalar() -> bool { true }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyType();
impl Coeff for PolyType {
    fn is_scalar() -> bool { false }
}


pub type SymbType = Option<&'static str>;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Poly<T: ScalarRing> {
    pub symb        : SymbType, // A literal for the indeterminates
    pub lead_scalar : T,
    pub terms : Vec<Monomial<T>>,
    // pub p_terms : Vec<Monomial<Poly<T>>>,
}

#[derive(Debug, Eq, PartialEq, Clone)]
enum CoeffType<T: ScalarRing> {
    Scalar(T),
    Poly(Poly<T>),
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Monomial<T: ScalarRing> {
    pub coeff : CoeffType<T>,
    pub deg   : usize,
}
// pub struct Poly<T: Ring, U: Coeff> {
//     pub symb  : SymbType, // A literal for the indeterminates
//     pub lead_scalar : <T>::BaseRing,
//     pub terms : Vec<Monomial<T, U>>,
// }

// Common aliases
// pub type PolyU<T> = Poly<T, ScalarType>;
// pub type PolyM<T> = Poly<T, PolyType>;

// #[derive(Debug, Eq, PartialEq, Clone)]
// pub struct Monomial<T: Ring, U: Coeff> {
//     pub coeff : T,
//     pub deg   : usize,
//     coeff_type  : PhantomData<U>,
// }

// Common aliases
// pub type MonomialU<T> = Monomial<T, ScalarType>;
// pub type MonomialM<T> = Monomial<T, PolyType>;

// I'm not planning on going more abstract than groups
impl<T: ScalarRing> Monomial<T> {
    pub fn new(coeff: CoeffType<T>, deg: usize) -> Monomial<T> {
        Monomial { coeff, deg }
    }
    pub fn zero() -> Self {
        Monomial::new(CoeffType::Scalar(<T>::zero()), 0 )
    } 
    pub fn one() -> Self {
        Monomial::new(CoeffType::Scalar(<T>::one()), 0)
    } 
    pub fn neg(&self) -> Self {
        Monomial::new(self.coeff.neg(), self.deg)
    }
}


impl<T: Ring, U: Coeff> Poly<T, U> {

    pub fn from_monomials(symb: SymbType, terms: Vec<Monomial<T, U>>) -> Result<Poly<T, U>, PolyErr> {
        // Coefficient vector cannot be empty
        // TODO Probably need to go through and sort the list as well
        // and remove duplicates
        Ok (
            Poly { 
                symb,
                lead_scalar: <T>::BaseRing::one(),
                terms,
            }
        )
    }

    // Assumes the last monomial has the greatest degree
    // TODO this fails in the multivariate case, unless we have a grelex
    // monomial ordering
    pub fn deg(&self) -> usize {
        self.terms[self.terms.len()-1].deg
    }
}

pub fn expand<T: Ring, U: Coeff>(input: &[Monomial<T, U>], n: usize) -> Vec<T> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<T> = Vec::with_capacity(n);
    for mono in input.into_iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mono.deg, <T>::zero());
        result.push(mono.coeff.clone());
    }
    // Pad the rest
    result.resize(n, <T>::zero());
    result
}


use std::fmt;

impl<T> fmt::Display for Poly<T, PolyType> 
    where T: fmt::Display + Ring {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "{}", 0)
        }

        let symb = if let Some(s) = self.symb {s} else {"~"};

        let mut mono_iter = self.terms.iter();

        let mut acc = if self.terms[0].deg == 0 {
            format!("{}", mono_iter.next().unwrap().coeff)
        } else {
            String::new()
        };

        for num in mono_iter {
            let exp = match num.deg {
                0 => String::new(),
                1 => format!("{}", symb),
                d => format!("{}^{}", symb, d),
            };
            acc.push_str(&format!(" + ({}){}", num.coeff, exp));
        }

        write!(f, "{}", acc)
    }

}

impl<T> fmt::Display for Poly<T, ScalarType> 
    where T: fmt::Display + Ring {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "{}", 0)
        }

        let symb = if let Some(s) = self.symb {s} else {"~"};

        let mut mono_iter = self.terms.iter();

        let mut acc = if self.terms[0].deg == 0 {
            format!("{}", mono_iter.next().unwrap().coeff)
        } else {
            String::new()
        };

        for num in mono_iter {
            let exp = match num.deg {
                0 => String::new(),
                1 => format!("{}", symb),
                d => format!("{}^{}", symb, d),
            };
            acc.push_str(&format!(" + {}{}", num.coeff, exp));
        }

        write!(f, "{}", acc)
    }

}

// impl std::str::FromStr for PolyU<ZZ> {
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

//         Ok(PolyU::from_monomials(None, acc)?)
//     }
// }