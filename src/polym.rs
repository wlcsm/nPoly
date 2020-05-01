// Multivariate polynomial implementation
//

// use std::fmt::Debug;
// use crate::algebras::Field;
use crate::algebras::polyring::*;
use crate::algebras::*;
use std::cmp::Ordering;
use generic_array::*;
use std::fmt::Debug;

// <><><><><><><><><><> Constructors <><><><><><><><><><> //
// impl<F, M> PRDomain<F, M> 
//     where F: Field, M: Multivariate {
//     // For now we have a cap of 6 variables, and so the total deg gets log2(6) = 3
//     // more bits than the others.
//     pub fn multivar(vars: Vec<String>) -> Self {
//         unimplemented!()
//     }
// }

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct MultiIndex<N: VarNumber> {
    total: usize,
    indices: GenericArray<usize, N>,
}

impl<N: VarNumber> MultiIndex<N> {
    fn new(ind: GenericArray<usize, N>) -> Self {
        MultiIndex {
            total: ind.iter().sum(),
            indices: ind,
        }
    }
}

impl<N: VarNumber + Debug> MultiIndexTrait for MultiIndex<N> {
    type N = N;

    fn lex(a: &Self, b: &Self) -> Ordering {
        for (i, j) in izip!(&a.indices, &b.indices) {
            if i != j {
                return i.cmp(&j)
            }
        }
        Ordering::Equal
    }
    fn glex(a: &Self, b: &Self) -> Ordering {
        match a.total.cmp(&b.total) {
            Ordering::Equal => <MultiIndex<N>>::lex(a, b),
            ord             => ord,
        }
    }
    fn grlex(a: &Self, b: &Self) -> Ordering {
        for (i, j) in izip!(&a.indices, &b.indices).rev() {
            if i != j {
                return j.cmp(&i) // Note that the order here is swapped
            }
        }
        Ordering::Equal
    }
}

impl<N: VarNumber> Zero for MultiIndex<N> {
    fn zero() -> Self { MultiIndex::new(GenericArray::default()) }
}

// and we'll implement IntoIterator
impl<N: VarNumber> IntoIterator for MultiIndex<N> {
    type Item = usize;
    type IntoIter = GenericArrayIter<usize, N>;

    fn into_iter(self) -> Self::IntoIter {
        self.indices.into_iter()
    }
}

use std::iter::FromIterator;

// and we'll implement FromIterator
impl<N: VarNumber> FromIterator<usize> for MultiIndex<N> {
    fn from_iter<I: IntoIterator<Item=usize>>(iter: I) -> Self {
        let mut arr: GenericArray<usize, N> = GenericArray::default();
        for (el, a) in iter.into_iter().zip(arr.iter_mut()) {
            *a = el;
        }
        MultiIndex::new(arr)
    }
}

use std::marker::PhantomData;

// TODO; Make this a macro to do all of them
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GLex<M: MultiIndexTrait>(PhantomData<M>);

impl<M: MultiIndexTrait> MonomialOrdering<M> for GLex<M> {
    fn cmp(a: &M, b: &M) -> Ordering { <M>::glex(a, b) }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Lex();

impl<M: MultiIndexTrait> MonomialOrdering<M> for Lex {
    fn cmp(a: &M, b: &M) -> Ordering { <M>::lex(a, b) }
}

pub trait MultiIndexTrait: Variate {
    type N: VarNumber;

    fn lex(a: &Self, b: &Self) -> Ordering;
    fn glex(a: &Self, b: &Self) -> Ordering;
    fn grlex(a: &Self, b: &Self) -> Ordering;
}

use std::cmp::{min, max};

impl<N: VarNumber> Variate for MultiIndex<N> {
    type NumVar = N;

    fn tot_deg(&self) -> usize {
        self.total
    }
    fn add(&self, other: &Self) -> Self {
        let new_indices = self.indices.iter()
            .zip(other.indices.iter())
            .map(|(a, b)| a + b)
            .collect();

        MultiIndex::new(new_indices)
    }
    fn sub(&self, other: &Self) -> Option<Self> {
        let mut new_indices = GenericArray::default();
        for i in 0 .. N::to_usize() {
            if self.indices[i] >= other.indices[i] {
                new_indices[i] = self.indices[i] - other.indices[i]
            } else {
                return None
            }
        }
        Some(MultiIndex::new(new_indices))
    }
    fn divides(&self, other: &Self) -> Option<bool> {
        Some(self.indices.iter()
                 .zip(other.indices.iter())
                 .all(|(a, b)| a >= b))
    }
    fn gcd(&self, other: &Self) -> Self {
        MultiIndex::new(
            self.indices.iter()
                .zip(other.indices.iter())
                .map(|(a, b)| *min(a, b)).collect())
    }
    fn lcm(&self, other: &Self) -> Self {
        MultiIndex::new(
            self.indices.iter()
                .zip(other.indices.iter())
                .map(|(a, b)| *max(a, b)).collect())
    }
}

// use std::fmt::*;
// use std::fmt;

// impl<> Debug for Multivariate<M> {
//     fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
//         write!(fmt, "")
//     }
// }