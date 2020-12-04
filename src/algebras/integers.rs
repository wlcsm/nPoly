extern crate alga;

use crate::algebras::*;
use crate::{impl_one, impl_zero};

use num_traits::identities::{One, Zero};

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct ZZ(pub i32);

impl ClosedAdd for ZZ {}
impl ClosedMul for ZZ {}
impl MyRing for ZZ {}


use std::ops;

impl_zero![ZZ, i32];
impl_one![ZZ, i32];

// Implement all the standard operations
impl_op_ex!(-   |a: &ZZ|         -> ZZ { ZZ(-a.0) });

impl_op_ex!(+   |a: &ZZ, b: &ZZ| -> ZZ { ZZ(a.0 + b.0) });
impl_op_ex!(-   |a: &ZZ, b: &ZZ| -> ZZ { ZZ(a.0 - b.0) });
impl_op_ex!(*   |a: &ZZ, b: &ZZ| -> ZZ { ZZ(a.0 * b.0) });

impl_op_ex!(+= |a: &mut ZZ, b: &ZZ|    { a.0 += b.0 });
impl_op_ex!(-= |a: &mut ZZ, b: &ZZ|    { a.0 -= b.0 });
impl_op_ex!(*= |a: &mut ZZ, b: &ZZ|    { a.0 *= b.0 });

impl ScalarRing for ZZ {
    const REGEX: &'static str = r"-?\d+";
}

use std::fmt;

impl fmt::Display for ZZ {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::str::FromStr for ZZ {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse::<i32>().map(|x| ZZ(x))
    }
}

impl EuclidDiv for ZZ {
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)> {
        if other.is_zero() {
            None
        } else {
            Some((
                ZZ(self.0.div_euclid(other.0)),
                ZZ(self.0.rem_euclid(other.0)),
            ))
        }
    }
}

impl EuclideanDomain for ZZ {
    fn divides(&self, other: &Self) -> Option<bool> {
        other
            .0
            .checked_rem_euclid(self.0)
            .and_then(|r| Some(r == 0))
    }
    fn gcd(&self, other: &Self) -> Self {
        if self.0.is_zero() {
            *other
        } else {
            ZZ(other.0 % self.0).gcd(&self)
        }
    }
    fn lcm(&self, other: &Self) -> Self {
        ZZ((self.0 * other.0) / self.gcd(&other).0)
    }
}
