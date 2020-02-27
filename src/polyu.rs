use crate::error::PolyErr;
use crate::algebras::*;

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
    pub fn zero() -> Self {
        Monomial { coeff: <T>::zero(), deg: 0 }
    } 
}

impl<T: Ring> Monomial<T> {
    pub fn one() -> Self {
        Monomial { coeff: <T>::one(), deg: 0 }
    } 
}

impl<T: Group> Monomial<T> {
    pub fn new(coeff: T, deg: usize) -> Monomial<T> {
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

    pub fn from_monomials(symb: Option<String>, terms: Vec<Monomial<T>>) -> Result<PolyU<T>, PolyErr> {
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
        self.terms[self.terms.len() - 1].deg
    }
}


impl<T: Group> AlgNeg for Monomial<T> {
    fn neg(&self) -> Self {
        Monomial::new(self.coeff.neg(), self.deg)
    }
}


// impl fmt::Display for PolyU<ZZ> {

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