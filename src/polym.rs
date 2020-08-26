/// Multivariate polynomial.
///
/// Most importantly it defines the MultiIndex struct which holds terms with multiple
/// indeterminates

use crate::algebras::polyring::*;
use crate::algebras::*;
use generic_array::*;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::marker::PhantomData;

/// Holds the indices of the term.
/// Also holds the total degree.
/// TODO Make it so that the first term in the monomial is the total degree.
#[derive(Clone, Debug, Hash)]
pub struct MultiIndex<N: VarNumber> {
    total: usize,
    indices: GenericArray<usize, N>,
}

impl<N:VarNumber> PartialEq for MultiIndex<N> {
    // First checks the total, then if it is true, checks the rest of the indices
    fn eq(&self, other: &Self) -> bool {
        self.total == other.total && self.indices == other.indices
    }
}
impl<N:VarNumber> Eq for MultiIndex<N> {}

/// Constructor
impl<N:VarNumber> MultiIndex<N> {
    fn new(indices: GenericArray<usize, N>) -> Self {
        MultiIndex {
            total: indices.iter().sum(),
            indices,
       }
    }
}

// TODO; Make this a macro to do all of them
#[derive(Clone, Eq, PartialEq, Debug, Hash)]
pub struct GLex<I: Monomial>(PhantomData<I>);

impl<I: Monomial> MonOrd for GLex<I> {
    type Index = I;

    fn cmp(itema: &I, itemb: &I) -> Ordering {
        match itema.tot_deg().cmp(&itemb.tot_deg()) {
            Ordering::Equal => itema.lex(itemb),
            ord => ord,
        }
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Lex();

use std::cmp::{max, min};

impl<N:VarNumber> ClosedAdd for MultiIndex<N> {}

pub fn binary_monomial_map<N: VarNumber>(lhs: &MultiIndex<N>, 
                                rhs: &MultiIndex<N>, 
                                map: fn(usize, usize) -> usize) -> MultiIndex<N> {
    MultiIndex::new(
        lhs.indices
            .iter()
            .zip(rhs.indices.iter())
            .map(|(a, b)| map(*a, *b))
            .collect()
    )
}

use std::ops::{Add, AddAssign};
impl<N:VarNumber> Add for MultiIndex<N> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        binary_monomial_map(&self, &other, |a, b| a + b)
    }
}

use num_traits::Zero;
impl<N:VarNumber> Zero for MultiIndex<N> {
    fn zero() -> Self {
        Self::new(GenericArray::default())
    }
    fn is_zero(&self) -> bool {
        self.total.is_zero()
    }
}

impl<N:VarNumber> AddAssign<Self> for MultiIndex<N> {
    fn add_assign(&mut self, other: Self) {
        for (a, b) in self.indices.iter_mut().zip(other.indices.iter()) {
            *a += b
        }
    }
}

impl<N:VarNumber> MyAddMonoid for MultiIndex<N> {
    fn ref_add(&self, other: &Self) -> Self {
        binary_monomial_map(self, other, |a, b| a + b)
    }
}


impl<N: VarNumber> Monomial for MultiIndex<N> {
    
    type NumVar = N;

    fn tot_deg(&self) -> usize {
        self.total
    }

    fn get(&self, ind: usize) -> Option<&usize> {
        self.indices.get(ind)
    }

    fn set(&mut self, ind: usize, val: usize) -> Option<()> {
        self.total += val - self.indices.get(ind)?;
        self.indices[ind] = val;
        Some(())
    }

    fn div(&self, other: &Self) -> Option<Self> {
        // FIXME returns None if it isn't divisible or if other is zero, these
        // two cases should be separate
        let mut res = MultiIndex::zero();
        let mut i = 0;
        for (a, b) in self.indices.iter().zip(other.indices.iter()) {
            if b <= a {
                res.indices[i] = a - b
            } else {
                return None;
            }
            i += 1;
        }
        Some(res)
    }

    /// The find map short-circuits when the two indices are not equal.
    /// If they are all equal then it returns None, in which case we 
    /// use a "map_or" to convert the None into a Ordering::Equal
    fn lex(&self, other: &Self) -> Ordering {
        izip!(&self.indices, &other.indices).find_map(|(a, b)| 
            match a.cmp(&b) {
                Ordering::Equal => None,
                ord             => Some(ord)
            })
            .map_or(Ordering::Equal , |ord| {ord})
    }

    fn gcd(&self, other: &Self) -> Self {
        binary_monomial_map(self, other, |a, b| min(a, b))
    }
    fn lcm(&self, other: &Self) -> Self {
        binary_monomial_map(self, other, |a, b| max(a, b))
    }
}

use std::fmt;
use crate::display::*;

// Problem is that it's hard to put an ordering on the coefficients because in finite fields
// thats quite ambiguous. I need it in the "if x < 0" line
// This will eventually have to be overcome some time.
impl<'a, P: PolyRing> fmt::Display for Poly<'a, P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Because we don't want a potential "+" out the front of the first term
        if self.is_zero() {
            write!(f, "{}", <P::Coeff>::zero())
        } else {
            let mut acc: String = show_term(&self.terms[0], &self.ring);

            self.terms
                .iter()
                .skip(1)
                .for_each(|x| acc.push_str(&format!(" + {}", show_term(x, &self.ring))));

            write!(f, "{}", acc)
        }
    }
}

mod kronecker {
    use generic_array::typenum::U1;
    use super::*;

    /// This will be used to hold the data to perform Kronecker substitution
    struct Kronecker_Data<N: VarNumber>(GenericArray<usize, N>);

    fn to_sparse_vec<P: PolyRing>(input: &DenseVec<P::Coeff>, kron_data: Option<Kronecker_Data<P::NumVar>>) -> Vec<Term<P>> {

        unimplemented!()
        // match kron_data {
        //     Some(_) => panic!("Can't do Kronecker substitution yet"),
        //     None    => {
        //         let terms = input.0
        //             .into_iter()
        //             .enumerate()
        //             .filter(|(_, c)| !c.is_zero())
        //             .map(|(i, c)| Term::new(c, UniIndex(i)))
        //             .collect();

        //         Poly::from_terms_unchecked(terms, Some(ring))
        //     }
        // }
    }
}
