use crate::algebras::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct ZZ(pub i32);
extern crate num_complex;

impl Ring for ZZ {
    type BaseRing = ZZ;
    fn is_poly() -> bool { false }

    fn add(&self, other: &Self) -> Self {
        ZZ(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        ZZ(self.0 - other.0)
    }
    fn neg(&self) -> Self {
        ZZ(-self.0)
    }
    fn zero() -> Self {
        ZZ(0)
    }
    fn mul(&self, other: &Self) -> Self {
        ZZ(self.0 * other.0)
    }
    fn one() -> Self {
        ZZ(1)
    }
}

impl ScalarRing for ZZ {
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