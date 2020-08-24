use crate::algebras::*;
use crate::{impl_zero, impl_one};

#[derive(Clone, Copy, Debug)]
pub struct RR(pub f64);

// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for RR {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).abs() < 0.000000001
    }
}

impl RR {
    pub fn new(n: f64) -> RR {
        RR(n)
    }
    /// TODO: Error checking?
    pub fn from_int(n: usize) -> RR {
        RR(n as f64)
    }
}


impl_zero![RR, f64];
impl_one![RR, f64];


use std::ops;

impl_op_ex!(-|a: &RR| -> RR { RR(-a.0) });

impl_op_ex!(+ |a: &RR, b: &RR| -> RR { RR(a.0 + b.0) });
impl_op_ex!(- |a: &RR, b: &RR| -> RR { RR(a.0 - b.0) });
impl_op_ex!(* |a: &RR, b: &RR| -> RR { RR(a.0 * b.0) });
impl_op_ex!(/ |a: &RR, b: &RR| -> RR { RR(a.0 / b.0) });

impl_op_ex!(+= |a: &mut RR, b: &RR| { a.0 += b.0 });
impl_op_ex!(-= |a: &mut RR, b: &RR| { a.0 -= b.0 });
impl_op_ex!(*= |a: &mut RR, b: &RR| { a.0 *= b.0 });
impl_op_ex!(/= |a: &mut RR, b: &RR| { a.0 /= b.0 });

impl ClosedAdd for RR {}
impl ClosedMul for RR {}

impl MyRing for RR {}
impl MyField for RR {}
impl ScalarField for RR {}

impl ScalarRing for RR {
    // Regex doesn't allow space on left or right hand side
    const REGEX: &'static str = r"-?\d*\.\d*\s*";
}

use std::fmt;

impl std::str::FromStr for RR {
    type Err = std::num::ParseFloatError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(RR(s.parse::<f64>()?))
    }
}

impl fmt::Display for RR {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
