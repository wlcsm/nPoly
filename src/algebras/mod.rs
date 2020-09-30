pub mod complex;
pub mod finite_field;
pub mod integers;
pub mod polyring;
pub mod real;
// use alga::general::{Ring, Field, Operator, Additive, Multiplicative};
use num_traits::{One, Zero};
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Alga package doesn't quite work for me because I need operations to
/// work when taking ownership or references

// The ref_add and ref_mul are additions and multiplications done through
// references
// The group trait is used in the MonomialIndex trait
pub trait ClosedAdd: Sized + Add<Self, Output = Self> {}

pub trait MyAddMonoid: Zero + Clone + ClosedAdd + PartialEq + Debug {
    fn ref_add(&self, other: &Self) -> Self;
}

impl<M: Zero + Sized + Copy + ClosedAdd + PartialEq + Debug> MyAddMonoid for M {
    fn ref_add(&self, other: &Self) -> Self {
        *self + *other
    }
}

pub trait ClosedMul: Sized + Mul<Self, Output = Self> {}

pub trait MyMulMonoid: One + Sized + Clone + ClosedMul + PartialEq + Debug {
    fn ref_mul(&self, other: &Self) -> Self;
}

impl<M: One + Copy + ClosedMul + PartialEq + Debug> MyMulMonoid for M {
    fn ref_mul(&self, other: &Self) -> Self {
        *self * *other
    }
}

pub trait MyAddGroup: MyAddMonoid + Neg<Output = Self> + Sub<Output = Self> {
    fn ref_sub(&self, other: &Self) -> Self;
}

impl<M: Copy + MyAddMonoid + Sub<Output = Self> + Neg<Output = Self>> MyAddGroup for M {
    fn ref_sub(&self, other: &Self) -> Self {
        *self - *other
    }
}

pub trait MyMulGroup: MyMulMonoid + Div<Output = Self> {
    fn ref_div(&self, other: &Self) -> Self;
}

impl<M: Copy + MyMulMonoid + Div<Output = Self>> MyMulGroup for M {
    fn ref_div(&self, other: &Self) -> Self {
        *self / *other
    }
}


#[macro_export]
macro_rules! derive_op {
    ($newtype:ty; $op:item $(,$gen_param:tt)?) => {
        impl$gen_param std::ops::$op for $newtype $gen_param {
            type Output = $newtype $gen_param;
            fn get_name!($op)(self, rhs: $alias) -> Self {
                self.0 get_op!($op) rhs.0
            }
        }
    };
}

macro_rules! impl_FooTrait {
    ($name:ty, $lifetime:tt) => {
        impl<$lifetime> $crate::FooTrait for $name<$lifetime> {  }
    };
    ($name:ty) => {
        impl $crate::FooTrait for $name {  }
    };
}

#[macro_export]
macro_rules! get_op {
    (Add) => {+};
    (Sub) => {-};
    (Mul) => {*};
}

#[macro_export]
macro_rules! get_name {
    (Add) => {add};
    (Sub) => {sub};
    (Mul) => {mul};
}

#[macro_export]
macro_rules! impl_zero {
    ($name:ident, $inner:ty) => {
        impl Zero for $name {
            #[inline]
            fn zero() -> Self {
                $name(<$inner>::zero())
            }
            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }
    };
}
#[macro_export]
macro_rules! impl_one {
    ($name:ident, $inner:ty) => {
        impl One for $name {
            #[inline]
            fn one() -> Self {
                $name(<$inner>::one())
            }
        }
    };
}

use std::ops::*;

// For the moment i'm also assuming that the rings are integral domains
pub trait MyRing: MyAddGroup + MyMulMonoid {}

use std::ops::{AddAssign, DivAssign, MulAssign};

pub trait ScalarRing:
    MyRing
    + Copy
    + std::fmt::Debug
    + std::str::FromStr
    + std::fmt::Display
    + MulAssign<Self>
    + AddAssign<Self>
    + SubAssign<Self>
{
    const REGEX: &'static str;
}

pub trait MyField: MyRing + MyMulGroup {}

pub trait ScalarField: MyField + ScalarRing + DivAssign<Self> {}

pub trait EuclidDiv: MyMulMonoid {
    /// In the form: self / other = Some(Quotient, Remainder)
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)>;
}

pub trait EuclideanDomain: MyRing {
    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    fn divides(&self, other: &Self) -> Option<bool>; // Option in case other == 0
}

impl<F: MyField> EuclideanDomain for F {
    #[inline(always)]
    fn gcd(&self, _other: &Self) -> Self {
        self.clone()
    }

    #[inline(always)]
    fn lcm(&self, _other: &Self) -> Self {
        self.clone()
    }

    #[inline(always)]
    fn divides(&self, _other: &Self) -> Option<bool> {
        if self.is_zero() {
            None
        } else {
            Some(true)
        }
    }
}

// /// A "qualifier" which imposes certain restrictions on the input
// /// At the moment I only have it for checking if something is nonzero
// /// But there are other things we could check here
// // pub struct Qual<T: MyAddMonoid>(T);

// impl<T: MyAddMonoid> Qual<T> {
//     fn new(item: T) -> Option<Qual<T>> {
//         if item.is_zero() {
//             None
//         } else {
//             Some(item)
//         }
//     }
// }
