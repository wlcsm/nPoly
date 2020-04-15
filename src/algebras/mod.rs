pub mod complex;
pub mod integers;
pub mod polyring;
pub mod finite_field;

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

pub trait ScalarRing: Ring + Copy + std::fmt::Debug {
    fn add_ass(&mut self, other: &Self);
    fn sub_ass(&mut self, other: &Self);
    fn mul_ass(&mut self, other: &Self);
}

pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}
// I have separated PolyMul from the PolyRing trait because the PolyRing trait
// can be implemented generically over univariate and multivariate polynomials.
// But the mul function cannot be (easily at least). This is because I need to implement
// Kronecker substitution or something. I might also just replace this trait with the 
// FastMul trait

// Ideally The PolyRing type should be the combination of these two, not PolyMul being the combination
// pub trait PolyMul: PolyRing {
//     fn mul(&self, other: &Self) -> Self;
// }

// pub trait PolyRing {
//     type BaseRing: ScalarRing + std::fmt::Debug;
//     fn zero(&self) -> Self;
//     fn is_zero(&self) -> bool;
//     fn add(&self, other: &Self) -> Self;
//     fn sub(&self, other: &Self) -> Self;
//     fn neg(&self) -> Self;
//     // Ring operations
//     fn is_one(&self) -> bool;
//     fn one(&self) -> Self;
//     fn scale(&self, scalar: Self::BaseRing) -> Self;
//     fn scale_ass(&mut self, scalar: Self::BaseRing);
// }