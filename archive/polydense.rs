use crate::algebras::*;
use std::ops::{Add, Sub, Mul, MulAssign};

#[derive(Debug, Clone, PartialEq)]
pub struct PolyDense<R: ScalarRing> {
    pub terms: Vec<R>,
}

impl<R: ScalarRing> PolyDense<R> {
    pub fn new(items: Vec<R>) -> PolyDense<R> {
        PolyDense { terms: items }
    }
    pub fn eval(&self, val: R) -> R {
        self.terms
            .iter()
            .rev()
            .fold(R::zero(), |acc, t| acc * val + *t)
    }
}

/// Dense Polynomial Arithmetic
///
/// Uses a macro to clean up the boilerplate
#[macro_export]
macro_rules! impl_arithmetic_polydense {
    ($op_trait:ident, $op_func:ident, $op_type:ty) => {
        impl<R: ScalarRing> $op_trait for $op_type {
            type Output = PolyDense<R>;

            fn $op_func(self, other: Self) -> PolyDense<R> {
                let new_coeff = self
                    .terms
                    .iter()
                    .zip(other.terms.iter())
                    .map(|(a, b)| a.$op_func(*b))
                    .collect();

                PolyDense { terms: new_coeff }
            }
        }
    };
}

impl_arithmetic_polydense![Add, add, PolyDense<R>];
impl_arithmetic_polydense![Add, add, &PolyDense<R>];
impl_arithmetic_polydense![Sub, sub, PolyDense<R>];
impl_arithmetic_polydense![Sub, sub, &PolyDense<R>];

/// Dense Polynomial Multiplication
///
/// Uses a simple implementation of the n^2 algorithm
impl<R: ScalarRing> Mul for PolyDense<R> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let d_a = self.terms.len() + 1;
        let d_b = rhs.terms.len() + 1;
        let d_res = d_a + d_b - 1;

        let mut res = vec![R::zero(); d_res];

        // We could try to reword this to do all the sums for a particular entry of "res" in one
        // go, rather than iterating over res multiply time. David Harvey does this in Sage
        for i in 0..d_a {
            for j in 0..d_b {
                res[i + j] += self.terms[i] * rhs.terms[j]
            }
        }

        // Remove any remaining zeros on the end
        while res.last().unwrap().is_zero() {
            res.pop();
        }

        PolyDense { terms: res }
    }
}

impl<R: ScalarRing> MulAssign for PolyDense<R> {
    fn mul_assign(&mut self, other: PolyDense<R>) {
        *self = self.clone() * other
    }
}

