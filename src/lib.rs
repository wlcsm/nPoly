/// The polynomial crate

// #[macro_use] extern crate log;

pub mod error;
pub mod polyu;
pub mod polym;
pub mod algebras;
pub mod fft;



#[cfg(test)]
mod tests {

    use crate::algebras::integers::ZZ;
    use crate::algebras::*;
    use super::polyu::*;

    #[test]
    fn basics() {
        let a = PolyU::from_coeff(None, vec![ZZ(1),ZZ(2),ZZ(3)]).unwrap();
        let b = PolyU::from_coeff(None, vec![ZZ(4),ZZ(5),ZZ(6)]).unwrap();
        let c = PolyU::from_coeff(None, vec![ZZ(0),ZZ(0),ZZ(1),ZZ(2)]).unwrap();

        // General adding with different lengths
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(5),ZZ(7),ZZ(9)]).unwrap(), a.add(&b));
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(4),ZZ(5),ZZ(7),ZZ(2)]).unwrap(), b.add(&c));

        // Testing the cleaning feature
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(0)]).unwrap(), b.sub(&b));

        // Negative numbers
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(-3),ZZ(-3),ZZ(-3)]).unwrap(), a.sub(&b));

        // Scaling
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(2),ZZ(4),ZZ(6)]).unwrap(), a.scale(ZZ(2)));
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(-2),ZZ(-4),ZZ(-6)]).unwrap(), a.scale(ZZ(-2)));
        assert_eq!(PolyU::from_coeff(None, vec![ZZ(0)]).unwrap(), a.scale(ZZ(0)));
    }
}
