/// The polynomial crate

extern crate generic_array;

pub mod error;
pub mod polyu;
pub mod polym;
pub mod algebras;
pub mod fast_mult;
pub mod fft;
mod mathutils;
pub mod sparse;

#[cfg(test)]
mod tests {

    use crate::algebras::integers::ZZ;
    use crate::algebras::*;
    use super::polyu::*;

    #[test]
    fn basics() {
        // TODO: Reimplement these tests. 

        // let a = Poly::<ZZ, Univariate>::from_coeff(None, vec![ZZ(1),ZZ(2),ZZ(3)]);
        // let b = Poly::from_coeff(None, vec![ZZ(4),ZZ(5),ZZ(6)]);
        // let c = Poly::from_coeff(None, vec![ZZ(0),ZZ(0),ZZ(1),ZZ(2)]);

        // // General adding with different lengths
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(5),ZZ(7),ZZ(9)]), a.add(&b));
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(4),ZZ(5),ZZ(7),ZZ(2)]), b.add(&c));

        // // Testing the cleaning feature
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(0)]), b.sub(&b));

        // // Negative numbers
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(-3),ZZ(-3),ZZ(-3)]), a.sub(&b));

        // Scaling
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(2),ZZ(4),ZZ(6)]).unwrap(), a.scale(ZZ(2)));
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(-2),ZZ(-4),ZZ(-6)]).unwrap(), a.scale(ZZ(-2)));
        // assert_eq!(Poly::from_coeff(None, vec![ZZ(0)]).unwrap(), a.scale(ZZ(0)));
    }
}
