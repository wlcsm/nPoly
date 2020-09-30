extern crate num_complex;

use crate::algebras::*;
use crate::{impl_one, impl_zero};
use num_complex::Complex64;
use std::ops;

#[derive(Clone, Copy, Debug)]
pub struct CC(pub Complex64);

impl CC {
    pub fn new(re: f64, im: f64) -> Self {
        Self(Complex64::new(re, im))
    }
}

// Approximations for now, I don't actually use this operation though in the FFT
impl PartialEq for CC {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).norm() < 0.000000001
    }
}

impl_zero![CC, Complex64];
impl_one![CC, Complex64];

// Implement all the standard operations
impl_op_ex!(-|a: &CC| -> CC { CC(-a.0) });

impl_op_ex!(+ |a: &CC, b: &CC| -> CC { CC(a.0 + b.0) });
impl_op_ex!(-|a: &CC, b: &CC| -> CC { CC(a.0 - b.0) });
impl_op_ex!(*|a: &CC, b: &CC| -> CC { CC(a.0 * b.0) });
impl_op_ex!(/ |a: &CC, b: &CC| -> CC { CC(a.0 / b.0) });

impl_op_ex!(+= |a: &mut CC, b: &CC| { a.0 += b.0 });
impl_op_ex!(-= |a: &mut CC, b: &CC| { a.0 -= b.0 });
impl_op_ex!(*= |a: &mut CC, b: &CC| { a.0 *= b.0 });
impl_op_ex!(/= |a: &mut CC, b: &CC| { a.0 /= b.0 });

impl ClosedAdd for CC {}
impl ClosedMul for CC {}

impl MyRing for CC {}
impl MyField for CC {}
impl ScalarField for CC {}

impl ScalarRing for CC {
    // Regex doesn't allow space on left or right hand side
    const REGEX: &'static str = r"-?\d*\.\d*\s*(\+|-)\s*\d*\.\d*i";
}

// // <><><><><> Constructors <><><><><> //
impl CC {
    pub fn from_re(val: i32) -> CC {
        CC(Complex64::new(val as f64, 0.0))
    }

    pub fn from_im(val: i32) -> CC {
        CC(Complex64::new(0.0, val as f64))
    }
}

use regex::Regex;

impl std::str::FromStr for CC {
    type Err = std::num::ParseFloatError;

    // TODO make the regex less ugly, also use the regex for RR in it
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let cc_regex =
            Regex::new(r"^(?P<re>-?\d*(?:\.\d*))\s*(?:\+|-)\s*(?:(?P<im>\d*(?:\.\d*))i)$").unwrap();
        assert!(cc_regex.is_match(s));
        let caps = cc_regex.captures(s).unwrap();
        Ok(CC(Complex64::new(
            caps["re"].parse::<f64>()?,
            caps["im"].parse::<f64>()?,
        )))
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn parsing_test() {
//         let a = "3.0 + 5.0i".parse::<CC>();
//         println!("{:?}", a);
//     }
// }

use std::fmt;

impl fmt::Display for CC {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

use crate::fft::SupportsFFT;
use std::f64::consts::PI;

// // <><><><><> Supports FFT Implementation <><><><><> //
impl SupportsFFT for CC {
    fn rou(n: usize, inv: bool) -> Vec<Self> {
        // Generates all the nth roots of unity
        // Changes it depending on whether computing the dft or the inverse
        let sign = if inv { 1.0 } else { -1.0 };
        let base = Complex64::new(0.0, sign * 2.0 * PI / n as f64);
        (0..n).map(|k| CC(base.scale(k as f64).exp())).collect()
    }

    fn divby2(&mut self, n: usize) {
        self.0 /= (1 << n) as f64
    }
}
