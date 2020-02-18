
extern crate num_complex;

use num_complex::Complex64;
use crate::algebras::*;


#[derive(Clone)]
pub struct CC(pub Complex64);

// Approximations for now, I don't acutally use this operation though in the FFT
impl PartialEq for CC {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).norm() < 0.000000001
    }
}

impl Eq for CC {}

impl Group for CC {
    const zero: Self = CC(Complex64::new(0.0, 0.0));
}

impl Mul for CC {
    type Output = Self;
     
    fn mul(self, other: Self) -> Self {
        CC(self.0 * other.0)
    }
}

impl MulAssign<Self> for CC {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl Div<i32> for CC {
    type Output = Self;

    fn div(self, rhs: i32) -> Self::Output {
        CC(self.0 / rhs as f64)
    }
}

impl DivAssign<i32> for CC {
    fn div_assign(self, rhs: i32) -> Self {
        CC(self.0 / rhs as f64)
    }
}

impl Neg for CC {
    type Output = CC;

    fn neg(self) -> Self::Output {
        CC(-self.0)
    }
}

impl Add for CC {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        CC(self.0 + other.0)
    }
}

impl Sub for CC {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        CC(self.0 - other.0)
    }
}

impl Ring for CC {
    const one: CC = CC(Complex64::new(1.0, 0.0));
}
