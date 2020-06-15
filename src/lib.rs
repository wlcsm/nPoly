#![feature(associated_type_bounds)]
#![feature(iterator_fold_self)]

// The associated type bound was used in the construction of the FPolyRing trait to 
// make it essentially a subclass of it
//
// The iterator_fold_self was used in the standard polynomial multiplication algorithm
// in which I used a fold_first function (this can probably be worked around)
/// The polynomial crate
extern crate generic_array;
#[macro_use]
extern crate itertools;
extern crate nalgebra;
extern crate num;

pub mod algebras;
pub mod error;
pub mod fast_mult;
pub mod fft;
pub mod ideals;
pub mod mathutils;
pub mod parse;
pub mod polym;
pub mod polyu;
pub mod sparse;

#[cfg(test)]
mod tests {

    use crate::algebras::integers::ZZ;
    use crate::algebras::polyring::*;
    use crate::algebras::*;

    #[test]
    fn basics() {
        let ring = PRDomain::new(vec!['x']);

        let a = Poly::from_coeff(&ring, vec![ZZ(1), ZZ(2), ZZ(3)]);
        let b = Poly::from_coeff(&ring, vec![ZZ(4), ZZ(5), ZZ(6)]);
        let c = Poly::from_coeff(&ring, vec![ZZ(0), ZZ(0), ZZ(1), ZZ(2)]);

        // General adding with different lengths
        assert_eq!(
            Poly::from_coeff(&ring, vec![ZZ(5), ZZ(7), ZZ(9)]),
            a.add(&b)
        );
        assert_eq!(
            Poly::from_coeff(&ring, vec![ZZ(4), ZZ(5), ZZ(7), ZZ(2)]),
            b.add(&c)
        );

        // Testing the cleaning feature
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(0)]), b.sub(&b));

        // Negative numbers
        assert_eq!(
            Poly::from_coeff(&ring, vec![ZZ(-3), ZZ(-3), ZZ(-3)]),
            a.sub(&b)
        );

        // Scaling
        assert_eq!(
            Poly::from_coeff(&ring, vec![ZZ(2), ZZ(4), ZZ(6)]),
            a.scale(ZZ(2))
        );
        assert_eq!(
            Poly::from_coeff(&ring, vec![ZZ(-2), ZZ(-4), ZZ(-6)]),
            a.scale(ZZ(-2))
        );
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(0)]), a.scale(ZZ(0)));
    }
}
