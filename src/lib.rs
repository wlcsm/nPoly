/// The polynomial structure
/// One core feature I'd like to enforce is that the vector in the polynomial
/// can never be empty

// #[macro_use] extern crate log;
#[macro_use] extern crate cached;

pub mod error;
pub mod polyu;
pub mod polym;
// use polym::PolyM;
use polyu::PolyU;
pub mod fft;

// #[derive(Debug)]
// pub enum Poly {
//     Multi(PolyM),
//     Uni(PolyU),
// }

// #[derive(Debug)]
// pub struct Multinomial {
//     norm : Monomial,
//     ext  : Poly,
// }


#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Monomial<T: Eq + PartialEq> {
    coeff : T,
    deg   : usize,
}

#[cfg(test)]
mod tests {

    use super::PolyU;
    // use super::fft::generate_rou;

    #[test]
    fn basics() {
        let a = PolyU::from_coeff("x".to_string(), vec![1,2,3]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![4,5,6]).unwrap();
        let c = PolyU::from_coeff("x".to_string(), vec![0,0,1,2]).unwrap();

        // General adding with different lengths
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![5,7,9]).unwrap(), a.add(&b));
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![4,5,7,2]).unwrap(), b.add(&c));

        // Testing the cleaning feature
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![0]).unwrap(), b.sub(&b));

        // Negative numbers
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![-3,-3,-3]).unwrap(), a.sub(&b));

        // Scaling
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![2,4,6]).unwrap(), a.scale(2));
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![-2,-4,-6]).unwrap(), a.scale(-2));
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![0]).unwrap(), a.scale(0));
    }

    #[test]
    fn multiplication() {
        let a = PolyU::from_coeff("x".to_string(), vec![1,1]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![1,2,1]).unwrap();

        // General multiplication
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![1,2,1]).unwrap(), a.mul(&a));
        assert_eq!(PolyU::from_coeff("x".to_string(), vec![1,3,3,1]).unwrap(), a.mul(&b));
    }
}
