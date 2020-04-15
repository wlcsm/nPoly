// Multivariate polynomial implementation
//

// use std::fmt::Debug;
// use crate::algebras::ScalarRing;
use crate::algebras::polyring::*;
use crate::algebras::*;
use std::cmp::Ordering;
use generic_array::*;
use std::marker::PhantomData;

// impl<T, M, N> PolyRing for PRDomain<T, MultiIndex>
//     where T: ScalarRing, M: MonomialOrdering, N: ArrayLength<usize> {

//     type Coeff = T;
//     type Var = Multivariate;
//     fn symb(&self) -> Vec<String> {
//         self.var.symb.into_iter().collect()
//     }
// }
// <><><><><><><><><><> Constructors <><><><><><><><><><> //
// impl<T, M> PRDomain<T, M> 
//     where T: ScalarRing, M: Multivariate {
//     // For now we have a cap of 6 variables, and so the total deg gets log2(6) = 3
//     // more bits than the others.
//     pub fn multivar(vars: Vec<String>) -> Self {
//         unimplemented!()
//     }
// }

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct MultiIndex<N: ArrayLength<usize>> {
    total: usize,
    indices: GenericArray<usize, N>,
}

impl<N: ArrayLength<usize> + Debug + Clone + Eq> MultiIndexTrait for MultiIndex<N> {
    fn lex(a: &Self, b: &Self) -> Ordering {
        unimplemented!()
    }
    fn glex(a: &Self, b: &Self) -> Ordering {
        a.total.cmp(&b.total)
    }
    fn grlex(a: &Self, b: &Self) -> Ordering {
        unimplemented!()
    }
}

impl<N: ArrayLength<usize>> Zero for MultiIndex<N> {
    fn zero() -> Self { MultiIndex { total: 0, indices: GenericArray::default() } }
}
impl<N: ArrayLength<usize> + Debug + Clone + Eq> IndexTrait for MultiIndex<N> {
    fn deg(&self) -> usize { self.total }
}

// TODO; Make this a macro to do all of them
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GLex();

impl MonomialOrdering for GLex {
    fn cmp(a: &Self::MInd, b: &Self::MInd) -> Ordering { <Self::MInd>::glex(a, b) }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Lex();

impl MonomialOrdering for Lex {
    fn cmp(a: &Self::MInd, b: &Self::MInd) -> Ordering { <Self::MInd>::lex(a, b) }
}

pub trait MonomialOrdering: Clone + Eq {
    type MInd: MultiIndexTrait;
    fn cmp(a: &Self::MInd, b: &Self::MInd) -> Ordering;
}

pub trait MultiIndexTrait: IndexTrait {
    type N: ArrayLength<usize>;

    fn lex(a: &Self, b: &Self) -> Ordering;
    fn glex(a: &Self, b: &Self) -> Ordering;
    fn grlex(a: &Self, b: &Self) -> Ordering;
    // Should put more arbitrary orderings here
}

#[derive(Clone, Eq, PartialEq)]
struct Multivariate<M: MonomialOrdering> {
    ord: PhantomData<M>,
}

impl<M> Variate for Multivariate<M> 
    where M: MonomialOrdering {

    type Index = M::MInd;

    fn cmp(a: &Self::Index, b: &Self::Index) -> Ordering {
        <M>::cmp(a, b)
    }
    fn tdeg(index: &Self::Index) -> usize {
        Self::Index::deg(index)
    }
    fn zero() -> Self::Index { Self::Index::zero() }
}

use std::fmt::*;
use std::fmt;

impl<M: MonomialOrdering> Debug for Multivariate<M> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        write!(fmt, "")
    }
}