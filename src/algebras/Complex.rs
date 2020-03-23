extern crate num_complex;

use num_complex::Complex64;
use crate::algebras::*;


#[derive(Clone, Copy, Debug)]
pub struct CC(pub Complex64);

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
// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for CC {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).norm() < 0.000000001
    }
}

impl Eq for CC {}

impl Ring for CC {
    type BaseRing = CC;
    fn is_poly() -> bool { false }

    fn add(&self, other: &Self) -> Self {
        CC(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        CC(self.0 - other.0)
    }
    fn neg(&self) -> Self {
        CC(-self.0)
    }
    fn zero() -> Self {
        CC(Complex64::new(0.0, 0.0))
    }
    fn mul(&self, other: &Self) -> Self {
        CC(self.0 * other.0)
    }
    fn one() -> Self {
        CC(Complex64::new(1.0, 0.0))
    }
}

impl ScalarRing for CC {
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