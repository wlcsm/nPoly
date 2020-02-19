extern crate num_complex;

use num_complex::Complex64;
use crate::algebras::*;


#[derive(Clone, Copy)]
pub struct CC(pub Complex64);

// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for CC {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).norm() < 0.000000001
    }
}

impl Eq for CC {}

impl AlgAdd for CC {
    fn add(&self, other: &Self) -> Self {
        CC(self.0 + other.0)
    }
}

impl AlgAddAssign for CC {
    fn add_ass(&mut self, other: &Self) {
        self.0 += other.0
    }
}

impl AlgSub for CC {
    fn sub(&self, other: &Self) -> Self {
        CC(self.0 - other.0)
    }
}

impl AlgSubAssign for CC {
    fn sub_ass(&mut self, other: &Self) {
        self.0 -= other.0
    }
}

impl AlgNeg for CC {
    fn neg(&self) -> Self {
        CC(-self.0)
    }
}

impl Group for CC {
    fn zero() -> Self {
        CC(Complex64::new(0.0, 0.0))
    }
}

impl AlgMul for CC {
    fn mul(&self, other: &Self) -> Self {
        CC(self.0 * other.0)
    }
}

impl AlgMulAssign for CC{
    fn mul_ass(&mut self, other: &Self) {
        self.0 *= other.0
    }
}

impl Ring for CC {
    fn one() -> Self {
        CC(Complex64::new(1.0, 0.0))
    }
    type BaseRing = Complex64;

    fn scale(&self, scalar: Self::BaseRing) -> Self {
        CC(self.0 * scalar)
    }

    fn scale_ass(&mut self, scalar: Self::BaseRing) {
        self.0 *= scalar
    }
}
