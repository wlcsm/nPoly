use crate::error::PolyErr;
use crate::algebras::*;

type SymbType = Option<String>;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyU<T: ScalarRing> {
    pub(crate) symb  : SymbType, // A literal for the indeterminates
    pub(crate) terms : Monomials<T>
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub(crate) struct Monomials<T: ScalarRing> {
    pub(crate) lead_scalar: T,
    coeffs : Vec<T>,
    degs   : Vec<usize>,
}

pub(crate) struct MonomialsIter<'a, T: ScalarRing> {
    data: &'a Monomials<T>,
    index: usize,
}

impl<'a, T: ScalarRing> Iterator for MonomialsIter<'a, T> {

    type Item = (T, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.index;
        match (self.data.coeffs.get(i), self.data.degs.get(i)) {
            (Some(c), Some(d)) => {
                self.index += 1;
                Some((*c, *d))
            },
            (_, _) => None
        }
    }
}
impl<'a, T: ScalarRing> Monomials<T> {
    pub(crate) fn iter(&'a self) -> MonomialsIter<'a, T> {
        MonomialsIter {
            data: self,
            index: 0
        }
    }
}


impl<T: ScalarRing> Monomials<T> {
    
    fn new((coeffs, degs): (Vec<T>, Vec<usize>)) -> Monomials<T> {
        // Should do more checks for duplicate elements and need to sort it as well
        if coeffs.len() != degs.len() {
            panic!("Attempting to make a Monomials struct with coefficients and \
                    degrees of different lengths")
        }

        Monomials {
            lead_scalar: <T>::one(),
            coeffs,
            degs
        }
    }
    pub fn with_capacity(n: usize) -> Monomials<T> {
        Monomials {
            lead_scalar : <T>::one(),
            coeffs      : Vec::with_capacity(n),
            degs        : Vec::with_capacity(n),
        }
    }
}

impl<T: ScalarRing> std::iter::FromIterator<(usize, T)> for Monomials<T> {

    fn from_iter<I: IntoIterator<Item=(usize, T)>>(iter: I) -> Self {
        // Note: It would be better if we had a 'with_capacity' call here
        let mut coeffs = Vec::new();
        let mut degs   = Vec::new();
        for (deg, coeff) in iter {
            coeffs.push(coeff);
            degs.push(deg);
        }
        Monomials::new((coeffs, degs))
    }
}

impl<T: ScalarRing> PolyU<T> {

    pub fn from_coeff(symb: SymbType, coeffs: Vec<T>) -> Result<PolyU<T>, PolyErr> {
        // Converts into a PolyU type. 
        // Automatically compress the terms argument

        let terms: Monomials<T> = coeffs.into_iter().enumerate()
                            .filter(|(_, c)| *c != <T>::zero())
                            .collect();

        Ok(PolyU { symb, terms })
    }

    pub(crate) fn from_terms(symb: SymbType, terms: Monomials<T>) -> Result<PolyU<T>, PolyErr> {
        // TODO Probably need to go through and sort the list as well and remove duplicates
        Ok (PolyU { symb, terms })
    }

    pub fn deg(&self) -> usize {
        // Covers the case of an empty vector as degree 0
        *self.terms.degs.last().unwrap_or(&0)
    }
}

impl<T: ScalarRing> Monomials<T> {

    // Only want these two arrays to be incremented at the same time
    pub fn push(&mut self, (coeff, deg): (T, usize)) {
        self.coeffs.push(coeff);
        self.degs.push(deg);
    }

    pub fn get(&self, i: usize) -> Option<(T, usize)> {
        if i < self.coeffs.len() {
            Some((self.coeffs[i], self.degs[i]))
        } else {
            None
        }
    }

    pub fn len(&self) -> usize { self.coeffs.len() }

    // A get_unchecked implementation
    pub fn get_uc(&self, i: usize) -> (T, usize) {
        unsafe{
            (*self.coeffs.get_unchecked(i), *self.degs.get_unchecked(i))
        }
    }
}

impl <T: ScalarRing> PolyU<T> {

    pub fn no_terms(&self) -> usize { self.terms.coeffs.len() }
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