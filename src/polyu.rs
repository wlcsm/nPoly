use std::fmt;
use crate::error::PolyErr;
use crate::Monomial;

#[derive(Debug, PartialEq)]
pub struct PolyU {
    pub symb  : String, // A literal for the indeterminates
    pub terms : Vec<Monomial>,
}

impl Monomial {
    fn new(coeff: i32, deg: usize) -> Result<Monomial, PolyErr> {
        // Need to impose some constraints
        Ok(Monomial { coeff, deg })
    }

    fn zero() -> Monomial {
        Monomial::new(0, 0).unwrap()
    }

    fn one() -> Monomial {
        Monomial::new(1, 0).unwrap()
    }
}

impl PolyU {

    pub fn from_one_vec_sparse(symb: String, terms: Vec<(i32, usize)>)
                               -> Result<PolyU, PolyErr> {

        // TODO Need to check some constraints
        // Coefficient vector cannot be empty
        match !terms.is_empty() {
            true  => Ok ( PolyU {
                    symb,
                    terms: terms.into_iter()
                                .map(|(a, b)| Monomial::new(a, b).unwrap())
                                .collect()
                    }),
            false => Err(PolyErr::EmptyPoly),
        }
    }

    pub fn from_coeff(symb: String, terms: Vec<i32>) -> Result<PolyU, PolyErr> {
        // Converts into a PolyU type. 
        // It does not accept empty vectors for the terms arguement.
        // It will automatically compress the terms argument
        // TODO Need to check some constraints
        // TODO Look into whether we should make some override not to compress
        //      for efficiency

        // Coefficient vector cannot be empty
        if terms.is_empty() {
            Err(PolyErr::EmptyPoly)
        } else {
            // Used in the filter_map method. If the coefficient is zero, it is 
            // filtered out, otherwise it is converted in to a monomial struct
            let to_monomial = |(i,c)| match c {
                0  => None,
                x  => Some(Monomial::new(x,i).unwrap()),
            };

            let mut result: Vec<Monomial> = terms.into_iter().enumerate()
                                                 .filter_map(to_monomial)
                                                 .collect();

            // Don't want empty vectors
            if result.is_empty() {
                result.push(Monomial::zero());
            }

            Ok(PolyU {
                symb,
                terms: result,
            })
        }
    }

    fn from(symb: String, terms: Vec<Monomial>) -> Result<PolyU, PolyErr> {

        // TODO Need to check some constraints
        // Coefficient vector cannot be empty
        match !terms.is_empty() {
            true  => Ok ( PolyU { symb, terms }),
            false => Err(PolyErr::EmptyPoly),
        }
    }

    // The following two functions should not require a symbol. Its a grey area
    pub fn zero(symb: String) -> PolyU {
        PolyU::from(symb, vec![Monomial::zero()]).unwrap()
    }

    pub fn one(symb: String) -> PolyU {
        PolyU::from(symb, vec![Monomial::one()]).unwrap()
    }

    pub fn deg(&self) -> usize {
        self.terms.iter().map(|x| x.deg).max().unwrap()
    }

    /// Unoptimised solution
    pub fn eval(&self, point: i32) -> i32 {
        self.terms.iter().map(|x| x.coeff * int_pow(point, x.deg)).sum()
    }
}

// A memoised implementation of calculating integer powers
cached!{
    INT_POW;
    // b: base, e: exp
    fn int_pow(b: i32, e: usize) -> i32 = {
        match e {
            0 => 1,
            1 => b,
            _ => int_pow(b, e >> 1) * int_pow(b, e >> 1) * int_pow(b, e % 2)
        }
    }
}

impl fmt::Display for PolyU {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sign = |x: i32| { if x < 0 {" - "} else {" + "} };
        // Formats a nomial: Assumes that num is not zero
        let nomial = |num: &Monomial| -> String {format!("{}{}",
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

impl std::str::FromStr for PolyU {
    /// The function to parse a string into a polynomial type
    type Err = PolyErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Clean and remove square brackets
        let poly_iter = s[1..].trim()
                         .trim_matches(|p| p == '[' || p == ']' )
                         .split(',');

        // Parse each element into i32.
        let mut acc: Vec<Monomial> = Vec::new();
        for (i, x) in poly_iter.enumerate() {
            acc.push(Monomial::new(x.parse::<i32>()?, i).unwrap())
        };

        Ok(PolyU::from(String::from("x"), acc)?)
    }
}

/// When benchmarking it would be interesting to see whether an
/// implementation with iterators (which would be much smaller code)
/// would be faster.
impl PolyU {
    // Generic because I think I might use this kind of operation a lot
    // Also allows easy generalisation to group elements
    fn elementwise_binary<F>(&self, other: &PolyU, func: F) -> PolyU
    where
        F: Fn(i32, i32) -> i32
    {
        let (smol, bigg) =
            if self.deg() > other.deg() {
                (&other.terms, &self.terms)
            } else {
                (&self.terms, &other.terms)
            };

        let mut result: Vec<Monomial> = Vec::with_capacity(bigg.len());

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
                        0 => (),
                        a => result.push(Monomial::new(a, x.deg).unwrap()),
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
            result.push(Monomial::zero())
        };

        PolyU::from(self.symb.clone(), result).unwrap()
    }

    pub fn add(&self, other: &PolyU) -> PolyU {
        PolyU::elementwise_binary(self, other, |a, b| a + b)
    }

    pub fn sub(&self, other: &PolyU) -> PolyU {
        PolyU::elementwise_binary(self, other, |a, b| a - b)
    }


    /// A map along the coefficients
    /// I havent done the thing which checks if zeros have been introduced
    fn elementwise_unary<F>(&self, func: F) -> PolyU
    where
        F: Fn(&i32) -> i32
    {
        let mut result: Vec<Monomial> = Vec::new();
        for el in self.terms.iter() {
            match func(&el.coeff) {
                0 => (),
                x => result.push(Monomial::new(x, el.deg).unwrap()),
            };
        };

        // Vector cannot be empty
        if result.len() == 0 {
            result.push(Monomial::zero())
        }

        PolyU::from(
            self.symb.clone(),
            result,
        ).unwrap()
    }

    pub fn scale(&self, scalar: i32) -> PolyU {
        self.elementwise_unary(|x| scalar * x)
    }

    pub fn expand(&self) -> PolyU {

        let mut result: Vec<Monomial> = Vec::new();

        let mut i: usize = 0;
        for el in self.terms.iter() {
            if el.deg > i {
                for j in i..el.deg {
                    result.push(Monomial::new(0, j).unwrap());
                };
                i = el.deg;
            } else {
                result.push(el.clone());
                i += 1;
            }
        };
        PolyU::from(self.symb.clone(), result).unwrap()
    }

    pub fn compress(&self) -> PolyU {
        PolyU::from(
            self.symb.clone(),
            self.terms.clone().into_iter().filter(|x| x.coeff != 0).collect()
        ).unwrap()
    }

    fn to_monomials(acc: &Vec<i32>) -> Vec<Monomial> {

        // Remove the extra zeros when converting to monomial vector
        let mut result: Vec<Monomial> = Vec::new();
        for i in 0..acc.len() {
            match acc[i] {
                0 => (),
                x => result.push(Monomial::new(x, i).unwrap()),
            };
        }
        if result.is_empty() {
            result.push(Monomial::zero())
        }
        result
    }

    // Standard O(n^2) implementation
    pub fn mul(&self, other: &PolyU) -> PolyU {
        let mut acc: Vec<i32> = vec![0; self.deg() + other.deg() + 1];

        let self_exp = &self.expand();
        let other_exp = &other.expand();

        for i in 0..=self.deg() {
            for j in 0..=other.deg() {
                acc[i + j] += self_exp.terms[i].coeff * other_exp.terms[j].coeff;
            }
        };


        PolyU::from(
            self.symb.clone(),
            PolyU::to_monomials(&acc)
        ).unwrap()
    }
}
