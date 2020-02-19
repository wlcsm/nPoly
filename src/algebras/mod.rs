pub mod complex;
pub mod integers;
pub mod polyring;

use std::ops::*;

/// The original motivation for making my own copies of these traits is so that 
/// I can borrow self rather than taking ownership

pub trait AlgAdd {
    fn add(&self, other: &Self) -> Self;
}

pub trait AlgAddAssign {
    fn add_ass(&mut self, other: &Self);
}

pub trait AlgSub {
    fn sub(&self, other: &Self) -> Self;
}

pub trait AlgSubAssign {
    fn sub_ass(&mut self, other: &Self);
}

pub trait AlgNeg {
    fn neg(&self) -> Self;
}

pub trait AlgMul {
    fn mul(&self, other: &Self) -> Self;
}

pub trait AlgMulAssign {
    fn mul_ass(&mut self, other: &Self);
}

// Both groups and rings are assumed to be commutative
pub trait Group: AlgNeg 
                + AlgAdd + AlgAddAssign 
                + AlgSub + AlgSubAssign 
                + Sized + Eq + Clone {
    fn zero() -> Self;
 }

pub trait Ring: Group + AlgMul + AlgMulAssign {
    type BaseRing;
    fn scale(&self, scalar: Self::BaseRing) -> Self;
    fn scale_ass(&mut self, scalar: Self::BaseRing);
    fn one() -> Self;
}

pub trait Field: Ring + Div<Output = Self> {}