pub mod complex;

use std::ops::*;

// Both groups and rings are assumed to be commutative
pub trait Group: Add<Output = Self> + Neg<Output = Self>
                 + Sub<Output = Self> + Sized + Eq + Clone {
    fn zero() -> Self;

 }

pub trait Ring: Group + Mul<Output = Self> + Mul<i32, Output=Self> + MulAssign<Self> + 
                MulAssign<i32> + Div<i32, Output=Self> + DivAssign<i32> {
    fn one() -> Self;
}

pub trait Field: Ring + Div<Output = Self> {}


impl Group for i32 {
    fn zero() -> i32 {0}
}

impl Ring for i32 {
    fn one() -> i32 {1}
}