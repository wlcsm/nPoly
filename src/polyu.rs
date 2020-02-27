use crate::error::PolyErr;
use crate::algebras::*;

type SymbType = Option<String>;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyU<T: ScalarRing> {
    pub symb  : SymbType, // A literal for the indeterminates
    pub lead_scalar : T,
    pub terms : Vec<Monomial<T>>,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Monomial<T: ScalarRing> {
    pub coeff : T,
    pub deg   : usize,
}

impl<T: ScalarRing> Monomial<T> {
    pub fn zero() -> Self {
        Monomial { coeff: <T>::zero(), deg: 0 }
    } 
    pub fn one() -> Self {
        Monomial { coeff: <T>::one(), deg: 0 }
    } 
    pub fn new(coeff: T, deg: usize) -> Monomial<T> {
        Monomial { coeff, deg }
    }
    pub fn neg(&self) -> Self {
        Monomial::new(self.coeff.neg(), self.deg)
    }
}

// impl<T: ScalarRing> Coeff<T> {
//     // Adds a scalar and a polynomial quickly
//     fn s_plus_p(scalar: T, poly: PolyU<T>) -> PolyU<T> {
//         poly.add(&PolyU::from_coeff(None, vec![scalar]).unwrap())
//     }
//     // Checks if a poly could be compressed into a scalar. However I think this will
//     // have more performance loss than simply not checking
//     fn from(poly: PolyU<T>) -> Coeff<T> {
//         if poly.terms.len() <= 1 {
//             if poly.terms[0].deg == 0 {
//                 if let Coeff::S(n) = poly.terms[0].coeff {
//                     return Coeff::S(n.mul(&poly.lead_scalar))
//                 }
//             }
//         }
//         Coeff::P(poly)
//     }
//     pub fn add_ass(&mut self, other: &Self) {
//         match (self, other) {
//            (Coeff::S(s1), Coeff::S(s2)) => s1.add_ass(s2),
//            (Coeff::S(s1), Coeff::P(p2)) => {self = &mut Coeff::P(Coeff::s_plus_p(*s1,*p2));},
//            (Coeff::P(p1), Coeff::S(s2)) => {self = &mut Coeff::P(Coeff::s_plus_p(*s2,*p1))},
//            (Coeff::P(p1), Coeff::P(p2)) => {self = &mut Coeff::P(p1.add(p2))},
//         }
//     }
// }

// impl<T: ScalarRing> Ring for Coeff<T> {

//     type BaseRing = T;
//     // This is not techincally correct
//     fn is_poly() -> bool { true }

//     fn zero() -> Self {Coeff::S(<T>::zero())}
//     fn one() -> Self {Coeff::S(<T>::one())}

//     fn add(&self, other: &Self) -> Self {
//         match (self, other) {
//            (Coeff::S(s1), Coeff::S(s2)) => Coeff::S(s1.add(&s2)),
//            (Coeff::S(s1), Coeff::P(p2)) => Coeff::P(Coeff::s_plus_p(*s1, *p2)),
//            (Coeff::P(p1), Coeff::S(s2)) => Coeff::P(Coeff::s_plus_p(*s2, *p1)),
//            (Coeff::P(p1), Coeff::P(p2)) => Coeff::P(p1.add(&p2))
//         }
//     }

//     fn sub(&self, other: &Self) -> Self {
//         self.add(&other.neg())
//     }

//     fn neg(&self) -> Self {
//         match self {
//             Coeff::S(s) => Coeff::S(s.neg()),
//             Coeff::P(p) => Coeff::P(p.neg()),
//         }
//     }

//     fn mul(&self, other: &Self) -> Self {
//         match (self, other) {
//            (Coeff::S(s1), Coeff::S(s2)) => Coeff::S(s1.mul(&s2)),
//            (Coeff::S(s1), Coeff::P(p2)) => Coeff::P(p2.scale(*s1)),
//            (Coeff::P(p1), Coeff::S(s2)) => Coeff::P(p1.scale(*s2)),
//            (Coeff::P(p1), Coeff::P(p2)) => Coeff::P(p1.mul(&p2))
//         }
//     }

// }


impl<T: ScalarRing> PolyU<T> {

    pub fn from_coeff(symb: SymbType, coeffs: Vec<T>) -> Result<PolyU<T>, PolyErr> {
        // Converts into a PolyU type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument

        if coeffs.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            let terms = coeffs.into_iter().enumerate()
                              .filter(|(_, c)| *c != <T>::zero())
                              .map(   |(i, c)| Monomial::new(c, i))
                              .collect();

            Ok(PolyU::<T> { symb, lead_scalar: <T>::one(), terms })
        }
    }

    pub fn from_monomials(symb: SymbType, terms: Vec<Monomial<T>>) -> Result<PolyU<T>, PolyErr> {
        // TODO Probably need to go through and sort the list as well
        // and remove duplicates
        Ok (PolyU { symb, lead_scalar: <T>::one(), terms })
    }

    pub fn deg(&self) -> usize {
        if self.terms.is_empty() {
            0
        } else {
            self.terms[self.terms.len() - 1].deg
        }
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