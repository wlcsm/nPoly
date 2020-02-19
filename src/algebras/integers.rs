use crate::algebras::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct ZZ(pub i32);

impl AlgAdd for ZZ {
    fn add(&self, other: &Self) -> Self {
        ZZ(self.0 + other.0)
    }
}

impl AlgAddAssign for ZZ {
    fn add_ass(&mut self, other: &Self) {
        self.0 += other.0
    }
}

impl AlgSub for ZZ {
    fn sub(&self, other: &Self) -> Self {
        ZZ(self.0 - other.0)
    }
}

impl AlgSubAssign for ZZ {
    fn sub_ass(&mut self, other: &Self) {
        self.0 -= other.0
    }
}

impl AlgNeg for ZZ {
    fn neg(&self) -> Self {
        ZZ(-self.0)
    }
}

impl Group for ZZ {
    fn zero() -> Self {ZZ(0)}
}

impl AlgMul for ZZ {
    fn mul(&self, other: &Self) -> Self {
        ZZ(self.0 * other.0)
    }
}

impl AlgMulAssign for ZZ{
    fn mul_ass(&mut self, other: &Self) {
        self.0 *= other.0
    }
}

impl Ring for ZZ {
    fn one() -> Self {ZZ(1)}
    type BaseRing = i32;
    fn scale(&self, scalar: Self::BaseRing) -> Self {
        ZZ(self.0 * scalar)
    }
    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        self.0 *= scalar
    }
}

use std::fmt;

impl fmt::Display for ZZ {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
