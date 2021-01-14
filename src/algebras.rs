use alga::general::{ClosedAdd, ClosedMul, ClosedSub, ClosedDiv};
use num_traits::{One, Zero};
use std::fmt::Debug;
use std::fmt::Display;

/// Macro for creating and then automatically implementing my algebra traits
#[macro_export]
macro_rules! create_trait {
    ($trait_name:ident; $first_depen:ident, $($extra_depen:ident),*) => {
        pub trait $trait_name: $first_depen $(+ $extra_depen)* {}
        impl<T: $first_depen $(+ $extra_depen)*> $trait_name for T {}
    };
}

create_trait![AddMonoid  ; Zero, ClosedAdd];
create_trait![MulMonoid  ; One, ClosedMul];
create_trait![AddGroup   ; AddMonoid, ClosedSub];
create_trait![MulGroup   ; MulMonoid, ClosedSub];
create_trait![Ring       ; AddGroup, MulMonoid];
create_trait![ScalarRing ; Eq, Ring, Copy, Debug, Display];
create_trait![Field      ; Ring, MulGroup];
create_trait![ScalarField; ScalarRing, ClosedDiv];

pub trait EuclideanDomain: Sized {
    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    /// In the form: self / other = Some(Quotient, Remainder)
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)>;
}
