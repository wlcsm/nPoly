pub mod complex;
pub mod integers;
pub mod polyring;
// pub mod polyring_u;
// pub mod polyring_m;

/// The original motivation for making my own copies of these traits is so that 
/// I can borrow self rather than taking ownership

// For the moment i'm also assuming that the rings are integral domains
pub trait Ring: Sized + Eq + Clone {
    type BaseRing: ScalarRing + std::fmt::Debug;
    fn is_poly() -> bool;
    // Group operations
    fn zero() -> Self;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    // Ring operations
    fn one() -> Self;
    fn mul(&self, oter: &Self) -> Self;
}

pub trait ScalarRing: Ring + Copy + std::fmt::Debug {
    fn add_ass(&mut self, other: &Self);
    fn sub_ass(&mut self, other: &Self);
    fn mul_ass(&mut self, other: &Self);
}

pub trait PolyRing: Ring {
    // fn get_symb(&self) -> SymbType;
    // fn set_symb(&mut self, symb: SymbType); 
    fn scale(&self, scalar: Self::BaseRing) -> Self;
    fn scale_ass(&mut self, scalar: Self::BaseRing);
}