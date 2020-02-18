use std::fmt;
use crate::error::PolyErr;
use crate::algebras::*;
use std::ops::{Add, Neg, Mul, Div, Sub};

use crate::algebras::Complex::*;


#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyU<T: Group> {
    pub symb  : String, // A literal for the indeterminates
    pub terms : Vec<Monomial<T>>,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Monomial<T: Group> {
    pub coeff : T,
    pub deg   : usize,
}

// I'm not planning on going more abstract than groups
impl<T: Group> Monomial<T> {
    const zero: Monomial<T> = Monomial::new(<T>::zero, 0);
}

impl<T: Ring> Monomial<T> {
    const one: Monomial<T> = Monomial::new(<T>::one, 0);
}

impl<T: Group> Monomial<T> {
    fn new(coeff: T, deg: usize) -> Monomial<T> {
        // Need to impose some constraints
        Monomial { coeff, deg }
    }
}

impl<T: Group> PolyU<T> {

    pub fn from_one_vec_sparse(symb: String, terms: Vec<(T, usize)>)
                               -> Result<PolyU<T>, PolyErr> {

        // TODO Need to check some constraints
        // Coefficient vector cannot be empty
        if terms.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            Ok ( PolyU {
                symb,
                terms: terms.into_iter()
                            .map(|(c, d)| Monomial::new(c, d))
                            .collect()
                })
        }
    }

    pub fn from_coeff(symb: String, coeffs: Vec<T>) -> Result<PolyU<T>, PolyErr> {
        // Converts into a PolyU type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument

        // Coefficient vector cannot be empty
        if coeffs.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            // Used in the filter_map method. If the coefficient is zero, it is 
            // filtered out, otherwise it is converted in to a monomial struct

            let terms = Vec::new();
            for (i, c) in coeffs.into_iter().enumerate() {
                if c != <T>::zero {
                    terms.push(Monomial::new(c, i));
                }
            }
            // let mut terms: Vec<Monomial<T>> = coeffs.into_iter().enumerate().
            //                                         .filter(|(i, c)| c == 0)
            //                                      .filter_map()
            //                                      .collect();
            

            // Don't want empty vectors
            if terms.is_empty() {
                terms.push(Monomial::zero);
            }

            Ok(PolyU { symb, terms })
        }
    }

    fn from_monomials(symb: String, terms: Vec<Monomial<T>>) -> Result<PolyU<T>, PolyErr> {

        // TODO Need to check some constraints
        // Coefficient vector cannot be empty
        // TODO Probably need to go through and sort the list as well
        match !terms.is_empty() {
            true  => Ok ( PolyU::<T> { symb, terms }),
            false => Err(PolyErr::EmptyPoly),
        }
    }

    pub fn deg(&self) -> usize {
        self.terms.iter().map(|x| x.deg).max().unwrap()
    }

    // /// Unoptimised solution
    // pub fn eval(&self, point: i32) -> i32 {
    //     self.terms.iter().map(|x| x.coeff * int_pow(point, x.deg)).sum()
    // }
// // A memoised implementation of calculating integer powers
// cached!{
//     INT_POW;
//     // b: base, e: exp
//     fn int_pow(b: i32, e: usize) -> i32 = {
//         match e {
//             0 => 1,
//             1 => b,
//             _ => int_pow(b, e >> 1) * int_pow(b, e >> 1) * int_pow(b, e % 2)
//         }
//     }
// }

}

impl<T: Ring> Mul for &PolyU<T> {
    type Output = PolyU<T>;
     
    // Standard O(n^2) implementation
    fn mul(self, other: Self) -> PolyU<T> {
        let mut result = vec![<T>::zero; self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                result[a.deg + b.deg] = a.coeff + b.coeff;
            }
        }

        PolyU::from_coeff( self.symb.clone(), result ).unwrap()
    }
}

impl<T: Ring> Mul for PolyU<T> {
    type Output = Self;
     
    // Standard O(n^2) implementation
    fn mul(self, other: Self) -> Self {
        let mut result = vec![<T>::zero; self.deg() + other.deg() + 1];

        for a in self.terms.iter() {
            for b in other.terms.iter() {
                result[a.deg + b.deg] = a.coeff + b.coeff;
            }
        }

        PolyU::from_coeff( self.symb.clone(), result ).unwrap()
    }
}

impl<T: Ring> Mul<T> for PolyU<T> {
    type Output = Self;

    fn mul(self, scalar: T) -> PolyU<T> {

        let mut result: Vec<Monomial<T>> = Vec::new();
        for el in self.terms.iter() {
            match el.coeff * scalar {
                <T>::zero => (),
                x => result.push(Monomial::<T>::new(x, el.deg)),
            };
        };

        // Vector cannot be empty
        if result.len() == 0 {
            result.push(Monomial::<T>::zero)
        }

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
        let terms = self.terms.iter().map(|m| m.neg()).collect();
        PolyU::from_monomials("x".to_string(), terms).unwrap()
    }
}

impl<T: Group> Add for PolyU<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        PolyU::elementwise_binary(self, other, |a, b| a + b)
    }
}

impl<T: Group> Sub for PolyU<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        PolyU::elementwise_binary(self, other, |a, b| a - b)
    }
}

impl<T: Group> Group for PolyU<T> {
    // TODO the "x" here is a hack
    const zero: PolyU<T> =
        PolyU::from_monomials("x".to_string(), vec![Monomial::zero]).unwrap();
}

impl<T: Ring> Ring for PolyU<T> {
    // TODO the "x" here is a hack
    const one: PolyU<T> =
        PolyU::from_monomials("x".to_string(), vec![Monomial::one]).unwrap();
}

impl<T: Group> PolyU<T> {

    // Generic because I think I might use this kind of operation a lot
    // Also allows easy generalisation to group elements
    fn elementwise_binary<F>(self, other: Self, func: F) -> Self
    where
        F: Fn(T, T) -> T
    {
        let (smol, bigg) =
            if self.deg() > other.deg() {
                (&other.terms, &self.terms)
            } else {
                (&self.terms, &other.terms)
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
                    match func(x.coeff, y.coeff) {
                        <T>::zero => (),
                        a => result.push(Monomial::<T>::new(a, x.deg)),
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
            result.push(Monomial::<T>::zero)
        };

        PolyU::from_monomials(self.symb.clone(), result).unwrap()
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

        Ok(PolyU::<i32>::from_monomials(String::from("x"), acc)?)
    }
}

/// When benchmarking it would be interesting to see whether an
/// implementation with iterators (which would be much smaller code)
/// would be faster.
impl<T: Ring> PolyU<T> {

    // /// A map along the coefficients
    // /// I havent done the thing which checks if zeros have been introduced
    // fn elementwise_unary<F>(&self, func: F) -> PolyU<i32>
    // where
    //     F: Fn(&i32) -> i32
    // {
    //     let mut result: Vec<Monomial<i32>> = Vec::new();
    //     for el in self.terms.iter() {
    //         match func(&el.coeff) {
    //             0 => (),
    //             x => result.push(Monomial::<i32>::new(x, el.deg).unwrap()),
    //         };
    //     };

    //     // Vector cannot be empty
    //     if result.len() == 0 {
    //         result.push(Monomial::<i32>::zero())
    //     }

    //     PolyU::<i32>::from(
    //         self.symb.clone(),
    //         result,
    //     ).unwrap()
    // }


    pub fn expand(&self) -> PolyU<T> {

        let mut result: Vec<Monomial<T>> = Vec::new();

        let mut i: usize = 0;
        for el in self.terms.iter() {
            if el.deg > i {
                for j in i..el.deg {
                    result.push(Monomial::new(<T>::zero, j));
                };
                i = el.deg;
            } else {
                result.push(el.clone());
                i += 1;
            }
        };
        PolyU::from_monomials(self.symb.clone(), result).unwrap()
    }

    pub fn compress(&self) -> PolyU<T> {
        PolyU::from_monomials(
            self.symb.clone(),
            self.terms.clone().into_iter().filter(|x| x.coeff != <T>::zero).collect()
        ).unwrap()
    }
}
