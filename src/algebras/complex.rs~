extern crate num_complex;
extern crate alga;

use alga::general::{Multiplicative, Additive, Identity, AbstractMagma, TwoSidedInverse};
use crate::algebras::*;
// use crate::fft::SupportsFFT;
use num_complex::Complex64;
use std::f64::consts::PI;

#[derive(Alga)]
#[alga_traits(Field(Additive, Multiplicative))]
#[derive(Clone, Copy, Debug)]
pub struct CC(pub Complex64);

// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for CC {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).norm() < 0.000000001
    }
}

impl Identity<Additive> for CC {
    fn identity() -> CC {
        CC::from_re(0)
    }
}

impl AbstractMagma<Additive> for CC {
    fn operate(&self, rhs: &Self) -> CC {
        CC(self.0 + rhs.0)
    }
}
impl TwoSidedInverse<Additive> for CC {
    fn two_sided_inverse(&self) -> Self {
        CC(-self.0)
    }
}

impl Identity<Multiplicative> for CC {
    fn identity() -> CC {
        CC::from_re(1)
    }
}

impl AbstractMagma<Multiplicative> for CC {
    fn operate(&self, rhs: &Self) -> CC {
        CC(self.0 * rhs.0)
    }
}
impl TwoSidedInverse<Multiplicative> for CC {
    fn two_sided_inverse(&self) -> Self {
        CC(1.0 / self.0)
    }
}

impl std::cmp::Eq for CC {}

// // <><><><><> Constructors <><><><><> //
impl CC {
    pub fn new(val: Complex64) -> CC {
        CC(val)
    }
    pub fn from_re(val: i32) -> CC {
        CC(Complex64::new(val as f64, 0.0))
    }

    pub fn from_im(val: i32) -> CC {
        CC(Complex64::new(0.0, val as f64))
    }
}

// impl Zero for CC {
//     fn zero() -> Self {
//         CC(Complex64::new(0.0, 0.0))
//     }
// }

// impl One for CC {
//     fn one() -> Self {
//         CC(Complex64::new(1.0, 0.0))
//     }
// }

// // <><><><><> Ring Implementation <><><><><> //
// impl Ring for CC {
//     type BaseRing = CC;

//     fn add(&self, other: &Self) -> Self {
//         CC(self.0 + other.0)
//     }
//     fn sub(&self, other: &Self) -> Self {
//         CC(self.0 - other.0)
//     }
//     fn neg(&self) -> Self {
//         CC(-self.0)
//     }
//     fn mul(&self, other: &Self) -> Self {
//         CC(self.0 * other.0)
//     }
// }

// use regex::Regex;

// impl std::str::FromStr for CC {
//     type Err = std::num::ParseFloatError;

//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         let cc_regex =
//             Regex::new(r"^(?P<re>-?\d*(?:\.\d*))\s*(?:\+|-)\s*(?:(?P<im>\d*(?:\.\d*))i)$").unwrap();
//         assert!(cc_regex.is_match(s));
//         let caps = cc_regex.captures(s).unwrap();
//         Ok(CC(Complex64::new(
//             caps["re"].parse::<f64>()?,
//             caps["im"].parse::<f64>()?,
//         )))
//     }
// }
// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn parsing_test() {
//         let a = "3.0 + 5.0i".parse::<CC>();
//         println!("{:?}", a);
//     }
// }

// use std::fmt;

// impl fmt::Display for CC {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         write!(f, "{}", self.0)
//     }
// }

// // <><><><><> Scalar Ring Implementation <><><><><> //
// impl ScalarRing for CC {
//     // Regex doesn't allow space on left or right hand side
//     const REGEX: &'static str = r"-?\d*\.\d*\s*(\+|-)\s*\d*\.\d*i";
//     fn add_ass(&mut self, other: &Self) {
//         self.0 += other.0
//     }

//     fn sub_ass(&mut self, other: &Self) {
//         self.0 -= other.0
//     }
//     fn mul_ass(&mut self, other: &Self) {
//         self.0 *= other.0
//     }
// }

// impl EuclideanDomain for CC {
//     fn divides(&self, _other: &Self) -> Option<bool> {
//         if *self != CC::zero() {
//             Some(true)
//         } else {
//             None
//         }
//     }
//     fn gcd(&self, _other: &Self) -> Self {
//         if *self != CC::zero() {
//             CC::one()
//         } else {
//             CC::zero()
//         }
//     }
//     fn lcm(&self, other: &Self) -> Self {
//         CC((self.0 * other.0) / self.gcd(&other).0)
//     }
// }

// impl Field for CC {
//     fn div(&self, other: &Self) -> Option<CC> {
//         if *other == CC::zero() {
//             None
//         } else {
//             Some(CC(self.0 / other.0))
//         }
//     }
// }

// // <><><><><> Supports FFT Implementation <><><><><> //
// impl SupportsFFT for CC {
//     fn rou(n: usize, inv: bool) -> Vec<Self> {
//         // Generates all the nth roots of unity
//         // Changes it depending on whether computing the dft or the inverse
//         let sign = if inv { 1.0 } else { -1.0 };
//         let base = Complex64::new(0.0, sign * 2.0 * PI / n as f64);
//         (0..n).map(|k| CC(base.scale(k as f64).exp())).collect()
//     }

//     fn divby2(&mut self, n: usize) {
//         self.0 /= (1 << n) as f64
//     }
// }
