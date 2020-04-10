// Multivariate polynomial implementation
//

// use std::fmt::Debug;
// use crate::algebras::ScalarRing;
use crate::algebras::polyring::*;
use crate::algebras::*;
use std::cmp::Ordering;
use std::marker::PhantomData;

// <><><><><><><><><><> Constructors <><><><><><><><><><> //
impl<T: ScalarRing, M: MonomialOrdering> PRDomain<T, Multivariate<M>> {
    // For now we have a cap of 6 variables, and so the total deg gets log2(6) = 3
    // more bits than the others.
    pub fn multivar(vars: Vec<String>) -> Self {
        // Checks
        if vars.len() < 2 {
            panic!("Bro plz no")
        }

        // Generate Masks
        let n = vars.len();
        let per_var = 61 / (n + 1); // Number of bits each variable gets
        let var_mask = (1 << per_var + 1) - 1; // Makes "per_var" numbers of leading zeros

        // My elaborate way of making the mask for the total degree
        let tot_mask = Mask {
            bit_mask : (1 << 63) - 1 << (63 - per_var - 3),
            shift    : (63 - per_var - 3) as u64,
        };

        let var_masks = (1..n).map(|i| Mask {
                                            bit_mask : var_mask << i * per_var,
                                            shift    : (i * per_var) as u64,
                                        }).collect();

        let masks = Multivariate { tot_mask, var_masks, ord: PhantomData };

        PRDomain::new(masks)
    }
}

// The bit_mask is the mask itself, the shift tells how much you need to shift it 
// down (the shift isn't absolutely necessary its just a small optimisation)
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Mask {
    bit_mask : u64,
    shift    : u64
}

trait MonomialOrdering: Clone + Eq + PartialEq + std::fmt::Debug {
    fn cmp(vars: &Multivariate<Self>, a: &TermIndex, b: &TermIndex) -> Ordering;
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GLex();

impl MonomialOrdering for GLex {
    fn cmp(vars: &Multivariate<Self>, a: &TermIndex, b: &TermIndex) -> Ordering {
        match vars.tdeg(a).cmp(&vars.tdeg(b)) {
            Ordering::Equal => a.0.cmp(&b.0),
            ret             => ret,
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
struct Multivariate<M: MonomialOrdering> {
    tot_mask: Mask,
    var_masks: Vec<Mask>,
    ord: PhantomData<M>,
}

impl<M: MonomialOrdering> Multivariate<M> {
    fn remove_tot(&self, deg: &TermIndex) -> u64 {
        !self.tot_mask.bit_mask & deg.0 // Bit clear
    }
}

impl<M: MonomialOrdering> Variate for Multivariate<M> {

    fn cmp(&self, a: &TermIndex, b: &TermIndex) -> Ordering {
        <M>::cmp(&self, a, b)
    }
    fn tdeg(&self, index: &TermIndex) -> usize {
        ((self.tot_mask.bit_mask & index.0) >> self.tot_mask.shift) as usize
    }
    fn zero() -> TermIndex { TermIndex(0) }
}