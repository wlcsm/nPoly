pub mod complex;
pub mod finite_field;
pub mod integers;
pub mod polyring;
pub mod real;

/// The original motivation for making my own copies of the sub and add traits is so that
/// I can borrow self rather than taking ownership

// The group trait is used in the MonomialIndex trait
pub trait Group: Zero + Sized + Eq + Clone {
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
}

// For the moment i'm also assuming that the rings are integral domains
pub trait Ring: Zero + One + Sized + Eq + Clone {
    type BaseRing: ScalarRing + std::fmt::Debug;
    // Group operations
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    // Ring operations
    fn mul(&self, other: &Self) -> Self;
}

pub trait ScalarRing:
    Ring + Copy + std::fmt::Debug + std::str::FromStr + std::fmt::Display
{
    const REGEX: &'static str;
    fn add_ass(&mut self, other: &Self);
    fn sub_ass(&mut self, other: &Self);
    fn mul_ass(&mut self, other: &Self);
}

pub trait EuclideanDomain: Ring {
    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    fn divides(&self, other: &Self) -> Option<bool>; // Option in case other == 0
}

// TODO: Make a macro which implements EuclideanDomain for fields

pub trait Field: ScalarRing + EuclideanDomain {
    fn div(&self, other: &Self) -> Option<Self>;
}

pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}
