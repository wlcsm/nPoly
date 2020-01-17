/// The polynomial structure
/// One core feature I'd like to enforce is that the vector in the polynomial
/// can never be empty

use std::fmt;

#[derive(Debug, PartialEq)]
pub struct Poly(Vec<i32>);

impl Poly {
    // Cannot be empty
    pub fn new(data: Vec<i32>) -> Result<Poly, &'static str> {

        match !data.is_empty() {
            true  => Ok ( Poly(data.clone()).clean() ),
            false => Err("Vector to initialise polynomial cannot be empty"),
        }
    }

    pub fn zero() -> Poly {
        Poly(vec![0])
    }

    pub fn one() -> Poly {
        Poly(vec![1])
    }

    pub fn data(&self) -> &Vec<i32> {
        &self.0
    }

    pub fn deg(&self) -> usize {
        self.0.len() - 1
    }

    /// Very badly optimised solution
    pub fn eval(&self, point: i32) -> i32 {
        self.0.iter().enumerate().map(|(x, c)| c * int_power(point, x)).sum() 
    }
}

/// Very basic implementation of an exponentiation function.
/// Just because I couldn't easily find something else.
fn int_power(base: i32, exp: usize) -> i32 {
    match exp {
        0 => 1,
        1 => base,
        _ => {
            let a = int_power(base, exp >> 1);
            a * a * int_power(base, exp % 2)
        },
    }
}

impl fmt::Display for Poly {
    
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let sign = |x| { if x < 0 {"-"} else {"+"} };

        let mut acc = format!("{}", self.0[0]); // Guaranteed not to be empty
        if self.0.len() > 1 {
            acc.push_str(&format!(" {} {}x", sign(self.0[1]), self.0[1].abs()));
            for (i, el) in (&self.0).into_iter().enumerate().skip(2) {
                match el.abs() {
                    0 => continue,
                    1 => acc.push_str(&format!(" {} x^{}",   sign(*el),    i)),
                    x => acc.push_str(&format!(" {} {}x^{}", sign(*el), x, i)),
                }
            }
        }
        write!(f, "{}", acc)
    }
}


impl std::str::FromStr for Poly {
    /// The function to parse a string into a polynomial type
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> { 
        // Clean and remove square brackets
        let poly_iter = s[1..].trim()
                         .trim_matches(|p| p == '[' || p == ']' )
                         .split(',');

        // Parse each element into i32.
        let mut acc: Vec<i32> = Vec::new();
        for el in poly_iter {
            let a = el.parse::<i32>()
                     .expect(format!("Could not parse '{}; as i32", el).as_str());
            acc.push(a);
        }
        Ok(Poly::new(acc).expect("Vector cannot be empty").clean())
    }
}

/// Gets rid of any extra zeros at the end of the vector
/// Will always leave at least one element in the vector
impl Poly {

    pub fn is_empty(self) -> bool {
        self.0.len() == 0
    }

    fn clean(mut self) -> Self {
        let mut i = self.0.len() - 1;
        while self.0[i] == 0 && i > 0 {
            self.0.pop();
            i -= 1
        }
        self
    }
}

/// When benchmarking it would be interesting to see whether an 
/// implementation with iterators (which would be must smaller code) 
/// would be faster.
impl Poly {
    // Generic because I think I might use this kind of operation a lot
    // Also allows easy generialisation to group elements
    fn elementwise_binary<F>(&self, other: &Poly, func: F) -> Poly 
    where 
        F: Fn(i32, i32) -> i32
    {
        let (largest, max, min) = if self.deg() > other.deg() {
            (self, self.deg(), other.deg())
        } else {
            (other, other.deg(), self.deg())
        };

        let mut result: Vec<i32> = Vec::with_capacity(max+1);

        for i in 0..=min {
            result.push(func(self.0[i], other.0[i]));
        }
        for i in min+1..=max {
            result.push(largest.0[i]);
        };

        Poly::new(result).unwrap()
    }
    
    pub fn add(&self, other: &Poly) -> Poly {
        Poly::elementwise_binary(self, other, |a, b| a + b)
    }

    pub fn sub(&self, other: &Poly) -> Poly {
        Poly::elementwise_binary(self, other, |a, b| a - b)
    }
    

    /// A map along the coefficients
    /// Should ensure that the function is binary.
    fn elementwise_unary<F>(&self, func: F) -> Poly 
    where 
        F: Fn(&i32) -> i32
    {
        Poly(self.0.iter().map(func).collect())
    }

    pub fn scale(&self, scalar: i32) -> Poly {
        self.elementwise_unary(|x| scalar * x).clean()
    }

    pub fn mul(&self, other: &Poly) -> Poly {
        let mut result = vec![0; self.deg() + other.deg() + 1];

        println!("{}, {}", self.deg(), other.deg());
        println!("{}", result.len());
        for i in 0..=self.deg() {
            for j in 0..=other.deg() {
                result[i + j] += self.0[i] * other.0[j];
            } 
        };
        Poly::new(result).unwrap().clean()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let a = Poly::new(vec![1,2,3]).unwrap();
        let b = Poly::new(vec![4,5,6]).unwrap();
        let c = Poly::new(vec![0,0,1,2]).unwrap();

        // General adding with different lengths
        assert_eq!(Poly::new(vec![5,7,9]).unwrap(), a.add(&b));
        assert_eq!(Poly::new(vec![4,5,7,2]).unwrap(), b.add(&c));

        // Testing the cleaning feature
        assert_eq!(Poly::new(vec![0]).unwrap(), b.sub(&b));

        // Negative numbers
        assert_eq!(Poly::new(vec![-3,-3,-3]).unwrap(), a.sub(&b));

        // Scaling
        assert_eq!(Poly::new(vec![2,4,6]).unwrap(), a.scale(2));
        assert_eq!(Poly::new(vec![-2,-4,-6]).unwrap(), a.scale(-2));
        assert_eq!(Poly::new(vec![0]).unwrap(), a.scale(0));
    }

    #[test]
    fn multiplication() {
        let a = Poly::new(vec![1,1]).unwrap();
        let b = Poly::new(vec![1,2,1]).unwrap();

        // General multiplication
        assert_eq!(Poly::new(vec![1,2,1]).unwrap(), a.mul(&a));
        assert_eq!(Poly::new(vec![1,3,3,1]).unwrap(), a.mul(&b));

        // More zeros
        let c = Poly::new(vec![0,0,1]).unwrap();
        let d = Poly::new(vec![-1,1]).unwrap();
        assert_eq!(Poly::new(vec![0,0,1,1]).unwrap(), a.mul(&c));
        assert_eq!(Poly::new(vec![-1,0,1]).unwrap(), a.mul(&d));
    }
}




