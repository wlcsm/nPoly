use std::fmt;
use crate::error::PolyErr;
use crate::algebras::*;
use std::ops::{Add, Neg, Mul, Div, Sub, DivAssign, MulAssign};

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyU<T: Group> {
    pub symb  : Option<String>, // A literal for the indeterminates
    pub terms : Vec<Monomial<T>>,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Monomial<T: Group> {
    pub coeff : T,
    pub deg   : usize,
}

// I'm not planning on going more abstract than groups
impl<T: Group> Monomial<T> {
    fn zero() -> Self {
        Monomial { coeff: <T>::zero(), deg: 0 }
    } 
}

impl<T: Ring> Monomial<T> {
    fn one() -> Self {
        Monomial { coeff: <T>::one(), deg: 0 }
    } 
}

impl<T: Group> Monomial<T> {
    fn new(coeff: T, deg: usize) -> Monomial<T> {
        Monomial { coeff, deg }
    }
}

impl<T: Group> PolyU<T> {

    pub fn from_coeff(symb: Option<String>, coeffs: Vec<T>) -> Result<PolyU<T>, PolyErr> {
        // Converts into a PolyU type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument

        if coeffs.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            let mut terms = Vec::new();
            for (i, c) in coeffs.into_iter().enumerate() {
                if c != <T>::zero() {
                    terms.push(Monomial::new(c, i));
                }
            }

            // Don't want empty vectors
            if terms.is_empty() {
                terms.push(Monomial::zero());
            }

            Ok(PolyU { symb, terms })
        }
    }

    fn from_monomials(symb: Option<String>, terms: Vec<Monomial<T>>) -> Result<PolyU<T>, PolyErr> {
        // Coefficient vector cannot be empty
        // TODO Probably need to go through and sort the list as well
        // and remove duplicates
        if terms.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            Ok (PolyU { symb, terms })
        }
    }

    // Assumes the last monomial has the greatest degree
    // TODO this fails in the multivariate case, unless we have a grelex
    // monomial ordering
    pub fn deg(&self) -> usize {
        self.terms[self.terms.len()].deg
    }
}


impl<T: Ring> Mul for PolyU<T> {
    type Output = Self;
     
    // Standard O(n^2) implementation
    fn mul(self, other: Self) -> Self {
        let mut result = vec![<T>::zero(); self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                result[a.deg + b.deg] = result[a.deg + b.deg].clone() + a.coeff.clone() * b.coeff.clone();
            }
        }

        PolyU::from_coeff( self.symb.clone(), result ).unwrap()
    }
}

impl<T: Ring> Mul<i32> for PolyU<T> {
    type Output = Self;

    fn mul(self, scalar: i32) -> PolyU<T> {

        if scalar == 0 {
            return PolyU::zero()
        }

        let result = self.terms.into_iter()
                         .map(|x|
                            Monomial::new(x.coeff * scalar, x.deg)
                        ).collect();

        PolyU::<T>::from_monomials(
            self.symb.clone(),
            result,
        ).unwrap()
    }
}


impl<T: Group> Neg for Monomial<T> {
    type Output = Monomial<T>;

    fn neg(self) -> Self::Output {
        Monomial::new(-self.coeff, self.deg)
    }
}

impl<T: Group> Neg for PolyU<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let terms = self.terms.into_iter().map(|m| m.neg()).collect();
        PolyU::from_monomials(self.symb.clone(), terms).unwrap()
    }
}


impl<T: Group> Add for PolyU<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        PolyU::elementwise_binary(self, other, |a, b| a + b)
    }
}

impl<T: Group> Sub for PolyU<T> {
    type Output = PolyU<T>;

    fn sub(self, other: Self) -> Self::Output {
        PolyU::elementwise_binary(self, other, |a, b| a - b)
    }
}

impl<T: Group> Group for PolyU<T> {
    fn zero() -> Self {
        PolyU { symb: None, terms: vec![Monomial::zero()] }
    }
}

impl<T: Ring> DivAssign<i32> for PolyU<T> {

    fn div_assign(&mut self, scalar: i32) {
        for x in self.terms.iter_mut() {
            x.coeff /= scalar;
        }
    }
}

impl<T: Ring> Div<i32> for PolyU<T> {
    type Output = Self;

    fn div(self, scalar: i32) -> Self::Output {
        let new_terms = self.terms.into_iter()
                            .map(|x|
                                Monomial::new(x.coeff / scalar, x.deg)
                            ).collect();
        PolyU::from_monomials(self.symb.clone(), new_terms).unwrap()
    }
}

impl<T: Ring> MulAssign<i32> for PolyU<T> {

    fn mul_assign(&mut self, scalar: i32) {
        for x in self.terms.iter_mut() {
            x.coeff *= scalar;
        }
    }
}

impl<T: Ring> MulAssign<Self> for PolyU<T> {
    fn mul_assign(&mut self, other: Self) {
        let aux = self.clone();
        self.terms = (aux * other).terms;
    }
}

impl<T: Ring> Ring for PolyU<T> {
    fn one() -> Self {
        PolyU { symb: None, terms: vec![Monomial::one()] }
    }
}

impl<T: Group> PolyU<T> {
    // TODO This whole thing needs to be redone
    fn elementwise_binary<F>(polya: PolyU<T>, polyb: PolyU<T>, func: F) -> PolyU<T>
    where
        F: Fn(T, T) -> T
    {
        let (smol, bigg) =
            if polya.deg() > polyb.deg() {
                (&polyb.terms, &polya.terms)
            } else {
                (&polya.terms, &polyb.terms)
            };

        let mut result: Vec<Monomial<T>> = Vec::with_capacity(bigg.len());

        let mut i = 0;

        for el in smol.iter() {
            match (el, &bigg[i]) {
                (x, y) if x.deg <  y.deg => {result.push(x.clone())},
                (x, y) if x.deg >  y.deg => { while x.deg > bigg[i].deg {
                    i += 1;
                    result.push(bigg[i].clone())}
                },
                (x, y) if x.deg == y.deg => {
                    i += 1;
                    let a = func(x.coeff.clone(), y.coeff.clone());
                    if a == <T>::zero() {
                        result.push(Monomial::<T>::new(a, x.deg));
                    }
                },
                _ => unreachable!(),
            };
        }

        // Append any remaining terms to the result vector
        for j in i..bigg.len() {
            result.push(bigg[j].clone())
        }

        if result.len() == 0 {
            result.push(Monomial::<T>::zero())
        };

        PolyU::from_monomials(polya.symb.clone(), result).unwrap()
    }
}


impl fmt::Display for PolyU<i32> {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sign = |x: i32| { if x < 0 {" - "} else {" + "} };
        // Formats a nomial: Assumes that num is not zero
        let nomial = |num: &Monomial::<i32>| -> String {format!("{}{}",
            if num.coeff.abs() == 1 {"".to_string()}  else {format!("{}", num.coeff.abs())},
            if num.deg         == 1 {"x".to_string()} else {format!("x^{}", num.deg)}
        )};

        // Perform an extra check on the first element.
        // This is only one where degree can be zero.
        let mut acc: String =
            if self.terms[0].deg == 0 {
                format!("{}", self.terms[0].coeff)
            } else {
                nomial(&self.terms[0])
            };

        self.terms.iter().skip(1)
                         .for_each(|x|
                                   acc.push_str(&format!("{}{}", sign(x.coeff),
                                                         nomial(x))));

        write!(f, "{}", acc)
    }

}

impl std::str::FromStr for PolyU<i32> {
    /// The function to parse a string into a polynomial type
    type Err = PolyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Clean and remove square brackets
        let poly_iter = s[1..].trim()
                         .trim_matches(|p| p == '[' || p == ']' )
                         .split(',');

        // Parse each element into i32.
        let mut acc: Vec<Monomial<i32>> = Vec::new();
        for (i, x) in poly_iter.enumerate() {
            acc.push(Monomial::<i32>::new(x.parse::<i32>()?, i))
        };

        Ok(PolyU::<i32>::from_monomials(None, acc)?)
    }
}