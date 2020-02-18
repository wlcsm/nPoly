pub mod Complex;

use std::ops::{Add, Neg, Mul, Div, Sub, MulAssign, DivAssign};

// Both groups and rings are assumed to be commutative
// not sure how I would implement that, but it would probably be
// symbolic and rely on the algorithms, actually it probably wouldn't
// be a big deal
pub trait Group: Add<Output = Self> + Neg<Output = Self> 
                 + Sub<Output = Self> + Sized + Eq + Clone {
    const zero: Self;

 }

pub trait Ring: Group + Mul<Output = Self> + 
                MulAssign<Self> + Div<i32> + DivAssign<i32> {
    const one: Self;
}

pub trait Field: Ring + Div<Output = Self> {}


impl Group for i32 {
    const zero: i32 = 0;
}

impl Ring for i32 {
    const one: i32 = 1;
}