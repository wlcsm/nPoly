// Multivariate polynomial implementation
//
use crate::algebras::*;
use crate::polyu::*;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct PolyM<T: ScalarRing> {
    pub symb  : Option<String>, // A literal for the indeterminates
    pub lead_scalar : T,
    pub terms : Vec<Multinomial<T>>,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum MultiCoeff<T: ScalarRing> {
    U(PolyU<T>),
    M(PolyM<T>),
} 

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Multinomial<T: ScalarRing> {
    pub coeff : MultiCoeff<T>,
    pub deg   : usize,
}