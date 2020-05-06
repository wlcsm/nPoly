use crate::algebras::*;

#[derive(Clone, Copy, Debug)]
pub struct RR(pub f64);

// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for RR {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).abs() < 0.000000001
    }
}
impl std::cmp::Eq for RR {}

// <><><><><> Constructors <><><><><> //
impl RR {
    pub fn new(val: f64) -> RR {
        RR(val)
    }
}

impl Zero for RR {
    fn zero() -> Self {
        RR(0.0)
    }
}

impl One for RR {
    fn one() -> Self {
        RR(1.0)
    }
}

// <><><><><> Ring Implementation <><><><><> //
impl Ring for RR {
    type BaseRing = RR;

    fn add(&self, other: &Self) -> Self {
        RR(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        RR(self.0 - other.0)
    }
    fn neg(&self) -> Self {
        RR(-self.0)
    }
    fn mul(&self, other: &Self) -> Self {
        RR(self.0 * other.0)
    }
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

// <><><><><> Scalar Ring Implementation <><><><><> //
impl ScalarRing for RR {
    const REGEX: &'static str = r"\d*\.\d*";
    fn add_ass(&mut self, other: &Self) {
        self.0 += other.0
    }

    fn sub_ass(&mut self, other: &Self) {
        self.0 -= other.0
    }
    fn mul_ass(&mut self, other: &Self) {
        self.0 *= other.0
    }
}

impl EuclideanDomain for RR {
    fn divides(&self, _other: &Self) -> Option<bool> {
        if *self != RR::zero() {
            Some(true)
        } else {
            None
        }
    }
    fn gcd(&self, _other: &Self) -> Self {
        if *self != RR::zero() {
            RR::one()
        } else {
            RR::zero()
        }
    }
    fn lcm(&self, other: &Self) -> Self {
        RR((self.0 * other.0) / self.gcd(&other).0)
    }
}

impl Field for RR {
    fn div(&self, other: &Self) -> Option<RR> {
        if *other == RR::zero() {
            None
        } else {
            Some(RR(self.0 / other.0))
        }
    }
}
