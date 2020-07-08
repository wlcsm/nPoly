extern crate alga;

use crate::algebras::*;
use alga::general::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct ZZ(pub i32);

impl Identity<Additive> for ZZ {
    fn identity() -> Self {
        ZZ(0)
    }
}

impl AbstractMagma<Additive> for ZZ {
    fn operate(&self, rhs: &Self) -> ZZ {
        ZZ(self.0 + rhs.0)
    }
}
impl TwoSidedInverse<Additive> for ZZ {
    fn two_sided_inverse(&self) -> Self {
        ZZ(-self.0)
    }
}

impl Identity<Multiplicative> for ZZ {
    fn identity() -> ZZ {
        ZZ(1)
    }
}

impl AbstractMagma<Multiplicative> for ZZ {
    fn operate(&self, rhs: &Self) -> ZZ {
        ZZ(self.0 * rhs.0)
    }
}

impl AbstractSemigroup<Additive> for ZZ {}
impl AbstractMonoid<Additive> for ZZ {}
impl AbstractQuasigroup<Additive> for ZZ {}
impl AbstractLoop<Additive> for ZZ {}
impl AbstractGroup<Additive> for ZZ {}
impl AbstractGroupAbelian<Additive> for ZZ {}

impl AbstractSemigroup<Multiplicative> for ZZ {}
impl AbstractMonoid<Multiplicative> for ZZ {}

impl AbstractRing<Additive, Multiplicative> for ZZ {}
impl AbstractRingCommutative<Additive, Multiplicative> for ZZ {}


// impl One for ZZ {
//     fn one() -> Self {
//         ZZ(1)
//     }
// }
// impl Zero for ZZ {
//     fn zero() -> Self {
//         ZZ(0)
//     }
// }

// impl Ring for ZZ {
//     type BaseRing = ZZ;

//     fn add(&self, other: &Self) -> Self {
//         ZZ(self.0 + other.0)
//     }
//     fn sub(&self, other: &Self) -> Self {
//         ZZ(self.0 - other.0)
//     }
//     fn neg(&self) -> Self {
//         ZZ(-self.0)
//     }
//     fn mul(&self, other: &Self) -> Self {
//         ZZ(self.0 * other.0)
//     }
// }

use std::fmt;

impl fmt::Display for ZZ {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// impl ScalarRing for ZZ {
//     const REGEX: &'static str = r"-?\d+";

//     fn add_ass(&mut self, other: &Self) {
//         self.0 += other.0
//     }
//     fn sub_ass(&mut self, other: &Self) {
//         self.0 -= other.0
//     }
//     fn mul_ass(&mut self, other: &Self) {
//         self.0 *= other.0
//     }
// }

impl std::str::FromStr for ZZ {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.parse::<i32>() {
            Ok(n) => Ok(ZZ(n)),
            Err(e) => Err(e),
        }
    }
}

// impl EuclideanDomain for ZZ {
//     fn divides(&self, other: &Self) -> Option<bool> {
//         other
//             .0
//             .checked_rem_euclid(self.0)
//             .and_then(|r| Some(r == 0))
//     }
//     fn gcd(&self, other: &Self) -> Self {
//         if self.0 == 0 {
//             *other
//         } else {
//             ZZ(other.0 % self.0).gcd(&self)
//         }
//     }
//     fn lcm(&self, other: &Self) -> Self {
//         ZZ((self.0 * other.0) / self.gcd(&other).0)
//     }
// }
