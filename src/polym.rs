// Multivariate polynomial implementation
//

use crate::algebras::polyring::*;
use crate::algebras::*;
use generic_array::*;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Clone, Debug, Hash)]
pub struct MultiIndex<N: VarNumber> {
    total: usize,
    indices: GenericArray<usize, N>,
}

impl<N:VarNumber> PartialEq for MultiIndex<N> {
    // Checking the total isn't necessary, but it is an optimisation
    // unless your monomials are homogeneous
    fn eq(&self, other: &Self) -> bool {
        if self.total != other.total {
            false
        } else {
            self.indices == other.indices
        }
    }
}
impl<N:VarNumber> Eq for MultiIndex<N> {}

// impl<N:VarNumber> Ord for MultiIndex<N> {
//     fn cmp(&self, other: &Self) -> Ordering {
//         match <M>::id() {
//             MonOrdType::Lex => lex(&self, &other),
//             MonOrdType::GLex => 
//                 match self.tot_deg().cmp(&other.tot_deg()) {
//                     Ordering::Equal => lex(self, other),
//                     ord => ord,
//                 }
//         }
//     }
//     // /// Iterates through 'iter' and short-circuits when 'pred' is false. It also returns that
//     // /// element in the Some() struct. Otherwise returns None
//     // fn s_c_iter<I: Iterator<Item = J>, J: Copy>(iter: I, pred: fn(J) -> bool) -> Option<J> {
//     //     for i in iter {
//     //         if pred(i) {
//     //             return Some(i);
//     //         }
//     //     }
//     //     None
//     // }
// }

// impl<N:VarNumber> MultiIndex<N> {
//     fn lex(arga: &Self, argb: &Self) -> Ordering {
//         for (a, b) in arga.iter().zip(argb.iter()) {
//             if a != b {
//                 return a.cmp(&b)
//             }
//         }
//         Ordering::Equal
//     }
// }

// impl<N: VarNumber> Ord for MultiIndex<GLex, N> {
//     fn cmp(&self, other: &Self) -> Ordering {
//         match self.total.cmp(&other.total) {
//             Ordering::Equal => <MultiIndex<GLex, N>>::lex(self, other),
//             ord => ord,
//         }
//     }
// }

// impl<N: VarNumber> Ord for MultiIndex<Lex, N> {
//     fn cmp(&self, other: &Self) -> Ordering {
//         <MultiIndex<Lex, N>>::lex(self, other)
//     }
// }


// impl<N: VarNumber> IntoIterator for GradIndex<N> {
//     type Item = usize;
//     type IntoIter = GenericArrayIter<usize, N>;

//     fn into_iter(self) -> Self::IntoIter {
//         self.indices.into_iter()
//     }
// }

// fn lex<N: VarNumber>(inda: &GradIndex<N>, indb: &GradIndex<N>) -> Ordering {
//     match s_c_iter(izip!(&inda.indices, &indb.indices), |(x, y)| x != y) {
//         Some((x, y)) => x.cmp(&y),
//         None => Ordering::Equal,
//     }
// }

// /// The index where we store the degree alongside the
// /// individual index values
// #[derive(PartialEq, Eq, Debug, Clone)]
// pub struct GradIndex<N: VarNumber> {
//     total: usize,
//     indices: GenericArray<usize, N>,
// }
// #[derive(PartialEq, Eq, Debug, Clone)]
// pub struct GradIndex<const NUMVAR: usize> {
//     total: usize,
//     indices: [usize; NUMVAR],
// }
// impl<N: VarNumber> GradIndex<N> {
//     fn new(indices: GenericArray<usize, N>) -> Self {
//         GradIndex {
//             total: indices.iter().sum(),
//             indices,
//         }
//     }
// }

// impl<N: VarNumber> ClosedAdd for GradIndex<N> {}
// impl<N: VarNumber> Add for GradIndex<N> {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         GradIndex::new(
//             self.indices
//                 .iter()
//                 .zip(other.indices.iter())
//                 .map(|(a, b)| a + b)
//                 .collect(),
//         )
//     }
// }
// impl<N: VarNumber> Zero for GradIndex<N> {
//     fn zero() -> Self {
//         GradIndex::new(GenericArray::default())
//     }
//     fn is_zero(&self) -> bool {
//         self.total.is_zero()
//     }
// }

// impl<N: VarNumber> MyAddMonoid for GradIndex<N> {
//     fn ref_add(&self, other: &Self) -> Self {
//         GradIndex::new(
//             self.indices
//                 .iter()
//                 .zip(other.indices.iter())
//                 .map(|(a, b)| a + b)
//                 .collect(),
//         )
//     }
// }


// impl<N: VarNumber> Indices for GradIndex<N> {
//     type NumVar = N;
//     const NUMVAR: usize = N::to_usize();

//     fn tot_deg(&self) -> usize {
//         self.total
//     }
    
//     fn lex(&self, other: &Self) -> Ordering {
//         match s_c_iter(izip!(&self.indices, &other.indices), |(x, y)| x != y) {
//             Some((x, y)) => x.cmp(&y),
//             None => Ordering::Equal,
//         }
//     }

//     fn div(&self, other: &Self) -> Option<Self> {
//         // FIXME returns None if it isn't divisible or if other is zero, these
//         // two cases should be separate
//         unimplemented!()
//         // let mut res = MultiIndex::zero();
//         // let mut i = 0;
//         // for (a, b) in self.indices.iter().zip(other.indices.iter()) {
//         //     if b <= a {
//         //         res.indices[i] = a - b
//         //     } else {
//         //         return None;
//         //     }
//         //     i += 1;
//         // }
//         // Some(res)
//     }

//     fn divides(&self, other: &Self) -> Option<bool> {
//         if self.is_zero() {
//             None
//         } else {
//             Some(self.div(other).is_some())
//         }
//     }
// }


/// Constructor
/// TODO The input really should be an 
impl<N:VarNumber> MultiIndex<N> {
    fn new(indices: GenericArray<usize, N>) -> Self {
        MultiIndex {
            total: indices.iter().sum(),
            indices,
       }
    }
}

// impl<N:VarNumber> MultiIndexTrait for MultiIndex<N, M> {
//     const NUMVAR = N;

//     fn lex(a: &Self, b: &Self) -> Ordering {
//         match s_c_iter(izip!(&a.indices, &b.indices), |(x, y)| x != y) {
//             Some((x, y)) => x.cmp(&y),
//             None => Ordering::Equal,
//         }
//     }
//     fn glex(a: &Self, b: &Self) -> Ordering {
//         match a.total.cmp(&b.total) {
//             Ordering::Equal => <MultiIndex<N>>::lex(a, b),
//             ord => ord,
//         }
//     }
//     fn grlex(a: &Self, b: &Self) -> Ordering {
//         match a.total.cmp(&b.total) {
//             Ordering::Equal => {
//                 match s_c_iter(izip!(&a.indices, &b.indices).rev(), |(x, y)| x != y) {
//                     Some((x, y)) => x.cmp(&y),
//                     None => Ordering::Equal,
//                 }
//             }
//             ord => ord,
//         }
//     }
// }


// and we'll implement IntoIterator
impl<N:VarNumber> IntoIterator for MultiIndex<N> {
    type Item = usize;
    type IntoIter = GenericArrayIter<usize, N>;

    fn into_iter(self) -> Self::IntoIter {
        self.indices.into_iter()
    }
}

// use std::iter::FromIterator;

// // and we'll implement FromIterator
// impl<N:VarNumber> FromIterator<usize> for MultiIndex<N> {
//     fn from_iter<I: IntoIterator<Item = usize>>(iter: I) -> Self {
//         let mut arr: GenericArray<usize, N> = GenericArray::default();
//         for (el, a) in iter.into_iter().zip(arr.iter_mut()) {
//             *a = el;
//         }
//         MultiIndex::new(arr)
//     }
// }

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

impl<N:VarNumber> Add for MultiIndex<N> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        MultiIndex::new(
            self.indices
                .iter()
                .zip(other.indices.iter())
                .map(|(a, b)| a + b)
                .collect(),
        )
    }
}

impl<N:VarNumber> MyAddMonoid for MultiIndex<N> {
    fn ref_add(&self, other: &Self) -> Self {
        MultiIndex::new(
            self.indices
                .iter()
                .zip(other.indices.iter())
                .map(|(a, b)| a + b)
                .collect(),
        )
    }
}

impl<N:VarNumber> Zero for MultiIndex<N> {
    fn zero() -> Self {
        Self::new(GenericArray::default())
    }

    fn is_zero(&self) -> bool {
        self.total.is_zero()
    }
}

use std::ops::{Add, AddAssign};

impl<N:VarNumber> AddAssign<Self> for MultiIndex<N> {
    fn add_assign(&mut self, other: Self) {
        for (a, b) in self.indices.iter_mut().zip(other.indices.iter()) {
            *a += b
        }
    }
}

use num_traits::Zero;

impl<N: VarNumber> Monomial for MultiIndex<N> {
    
    type NumVar = N;

    fn tot_deg(&self) -> usize {
        self.total
    }

    fn get(&self, ind: usize) -> Option<&usize> {
        self.indices.get(ind)
    }

    fn set(&mut self, ind: usize, val: usize) -> Option<()> {
        let change = val - self.indices.get(ind)?;
        self.indices[ind] = val;
        self.total += change;
        Some(())
    }

    // fn sub(&self, other: &Self) -> Option<Self> {
    //     let mut new_indices = GenericArray::default();
    //     for i in 0..N::to_usize() {
    //         if self.indices[i] >= other.indices[i] {
    //             new_indices[i] = self.indices[i] - other.indices[i]
    //         } else {
    //             return None;
    //         }
    //     }
    //     Some(MultiIndex::new(new_indices))
    // }
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

    // fn divides(&self, other: &Self) -> Option<bool> {
    //     // Evaluates self | other, "self divides other"
    //     Some(
    //         self.indices
    //             .iter()
    //             .zip(other.indices.iter())
    //             .all(|(a, b)| a <= b),
    //     )
    // }
    fn gcd(&self, other: &Self) -> Self {
        MultiIndex::new(
            self.indices
                .iter()
                .zip(other.indices.iter())
                .map(|(a, b)| *min(a, b))
                .collect(),
        )
    }
    fn lcm(&self, other: &Self) -> Self {
        MultiIndex::new(
            self.indices.iter()
                .zip(other.indices.iter())
                .map(|(a, b)| *max(a, b))
                .collect(),
        )
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
