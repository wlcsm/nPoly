/// The polynomial crate

// #[macro_use] extern crate log;

pub mod error;
pub mod polyu;
pub mod polym;
pub mod algebras;
pub mod fft;




#[cfg(test)]
mod tests {

    use super::polyu::*;

    #[test]
    fn basics() {
        let a = PolyU::from_coeff(None, vec![1,2,3]).unwrap();
        let b = PolyU::from_coeff(None, vec![4,5,6]).unwrap();
        let c = PolyU::from_coeff(None, vec![0,0,1,2]).unwrap();

        // General adding with different lengths
        assert_eq!(PolyU::from_coeff(None, vec![5,7,9]).unwrap(), a.clone() + b.clone());
        assert_eq!(PolyU::from_coeff(None, vec![4,5,7,2]).unwrap(), b.clone() + c.clone());

        // Testing the cleaning feature
        assert_eq!(PolyU::from_coeff(None, vec![0]).unwrap(), b.clone() - b.clone());

        // Negative numbers
        assert_eq!(PolyU::from_coeff(None, vec![-3,-3,-3]).unwrap(), a.clone() - b.clone());

        // Scaling
        assert_eq!(PolyU::from_coeff(None, vec![2,4,6]).unwrap(), a.clone() * 2);
        assert_eq!(PolyU::from_coeff(None, vec![-2,-4,-6]).unwrap(), a.clone() * -2);
        assert_eq!(PolyU::from_coeff(None, vec![0]).unwrap(), a.clone() * 0);
    }
}
