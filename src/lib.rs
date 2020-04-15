/// The polynomial crate

extern crate generic_array;

pub mod error;
pub mod polyu;
pub mod polym;
pub mod algebras;
pub mod fast_mult;
pub mod mathutils;
pub mod fft;
pub mod ideals;
pub mod sparse;

#[cfg(test)]
mod tests {

    use crate::algebras::integers::ZZ;
    use crate::algebras::*;
    use crate::algebras::polyring::*;
    use super::polyu::*;

    #[test]
    fn basics() {
        
        let ring = PRDomain::univar(1);

        let a = Poly::from_coeff(&ring, vec![ZZ(1),ZZ(2),ZZ(3)]);
        let b = Poly::from_coeff(&ring, vec![ZZ(4),ZZ(5),ZZ(6)]);
        let c = Poly::from_coeff(&ring, vec![ZZ(0),ZZ(0),ZZ(1),ZZ(2)]);

        // General adding with different lengths
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(5),ZZ(7),ZZ(9)]), a.add(&b));
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(4),ZZ(5),ZZ(7),ZZ(2)]), b.add(&c));

        // Testing the cleaning feature
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(0)]), b.sub(&b));

        // Negative numbers
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(-3),ZZ(-3),ZZ(-3)]), a.sub(&b));

        // Scaling
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(2),ZZ(4),ZZ(6)]), a.scale(ZZ(2)));
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(-2),ZZ(-4),ZZ(-6)]), a.scale(ZZ(-2)));
        assert_eq!(Poly::from_coeff(&ring, vec![ZZ(0)]), a.scale(ZZ(0)));
    }
}
