use crate::algebras::polyring::*;
use crate::ideals::*;

use crate::algebras::EuclideanDomain;

/// An iterator which returns an iterator of
/// all unordered pairs of the original iterator
// struct Unord<I: Iterator>{
//     first: I,
//     second: I,
//     curr_el: Option<I::Item>,
// }

// /// TODO also need to take into account whether we want the trace
// /// that is, i < j, or i <= j
// impl<I: Iterator + Clone> Unord<I> {
//     fn new(iter: I) -> Unord<I> {
//         let el = iter.next();
//         Unord { 
//             first: iter, 
//             second: iter.clone(),
//             curr_el: el,
//         }
//     }
// }

// // Then, we implement `Iterator` for our `Counter`:

// impl<I: Iterator> Iterator for Unord<I> {

//     type Item = (I::Item, I::Item);

//     fn next(&mut self) -> Option<Self::Item> {
//         match self.second.next() {
//             Some(b) => Some(self.curr_el.unwrap(), b)),
//             None => match curr_el {
//                 None => (None, None),
//                 Some(a) => {
//                     self.curr_el = self.first.next();
//                     self.sec = 
//                     (a, b)
//                 }
//             }
//         }
//     }
// }


impl<'a, P: FPolyRing> Ideal<'a, P> {
    pub fn is_groebner_basis(&self) -> bool {
        let n = self.gens.len();
        iproduct!((0..n), (0..n))
            .filter(|(i, j)| i < j)
            .map(|(i, j)| s_poly(&self.gens[i], &self.gens[j]).reduce(&self))
            .all(|x| x.is_zero())
    }

    /// Standard implementation of Buchberger's Algorithm
    pub fn bb_algorithm(&self) -> Self {
        let mut gb = self.clone();
        let mut n = gb.gens.len();
        let mut pairs = unord_pairs_int(n);

        // Iterates through all the pairs until there are no more pairs to check
        while !pairs.is_empty() {
            let mut new_pairs = Vec::new();

            for (i, j) in pairs {
                let r = s_poly(&gb.gens[i], &gb.gens[j]).reduce(&gb);
                if !r.is_zero() {
                    // Add the remainder to the ideal and update the pairs that we
                    // need to check
                    if let Some(_) = gb.add(r) {
                        n += 1;
                        new_pairs.append(&mut (0..n - 1).map(|k| (k, n-1)).collect());
                    }
                }
            }
            pairs = new_pairs;
        }
        gb
    }

    /// Implementation of the Improved Buchberger's algorithm found at the end of
    /// "Ideals, Varieties, and Algorithms" by Cox, Little, and O'Shea 4th edition.
    /// Not fully implemented, the trick from Proposition 3 hasn't yet been implemented as it is
    /// non-trivial
    pub fn bb_algorithm_impr(&self) -> Self {
        let mut gb = self.clone();
        let mut n = gb.gens.len();
        let mut unchecked_pairs = unord_pairs_int(n);

        // Iterates through all the pairs until there are no more pairs to check
        while !unchecked_pairs.is_empty() {
            let mut new_pairs = Vec::new();

            for (i, j) in unchecked_pairs {
                if gb.gens[i].lt().gcd(&gb.gens[j].lt()) != Term::one() {
                    let r = s_poly(&gb.gens[i], &gb.gens[j]).reduce(&gb);
                    if !r.is_zero() {
                        // Add the remainder to the ideal and update the pairs that we
                        // need to check
                        gb.add(r);
                        n += 1;
                        new_pairs.append(&mut (0..n - 1).map(|k| (k, n)).collect());
                    }
                }
            }
            unchecked_pairs = new_pairs;
        }
        gb
    }
}
fn unord_pairs_int(n: usize) -> Vec<(usize, usize)> {
    iproduct!((0..n), (0..n)).filter(|(i, j)| i >= j).collect()
}

fn s_poly<'a, P: FPolyRing>(lhs: &Poly<'a, P>, rhs: &Poly<'a, P>) -> Poly<'a, P> {
    let lcm = lhs.lt().lcm(&rhs.lt());
    let a_newlead = lcm.div(lhs.lt()).unwrap();
    let b_newlead = lcm.div(rhs.lt()).unwrap();
    lhs.term_scale(&a_newlead)
        .sub(&rhs.term_scale(&b_newlead))
}


impl<'a, P: FPolyRing> Poly<'a, P> {
    fn reduce(&self, g: &Ideal<'a, P>) -> Self {
        self.divpolys(&g.gens).1
    }

    pub fn divpolys(&self, f: &Vec<Self>) -> (Vec<Self>, Self) {
        let mut p = self.clone();
        let mut r = self.zero(); // This is a zero polynomial in the same ring as f
        let mut q = vec![self.zero(); f.len()];

        'outer: while !p.is_zero() {
            for i in 0..f.len() {
                if let Some(quo) = p.lt().div(f[i].lt()) {
                    p = p.sub(&f[i].term_scale(&quo));
                    q[i] = q[i].add(&Poly::from_terms(vec![quo], &q[i].ring));
                    continue 'outer;
                }
            }
            // Executes if no division occurred
            r = r.add(&Poly::from_terms(vec![p.terms.pop().unwrap()], &self.ring));
        }
        (q, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::RR;
    use crate::parse::*;
    use crate::polym::*;
    use crate::polyu::*;
    use generic_array::typenum::U2;

    /// Tests the Euclidean division algorithm for polynomials with multiple divisors
    #[test]
    fn division_test() {
        // Univariate
        let ring = PRDomain::<RR, UniIndex, UnivarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);

        // Multivariate
        let ring = PRDomain::<RR, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "3.0x^2y^6 + 5.0x^98y^2").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        println!("a = {}", a);
        println!("b = {}", b);
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);
    }

    /// Tests the normal Buchberger's Algorithm
    #[test]
    fn bb_alg_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]).bb_algorithm();
        assert!(r.is_groebner_basis());
    }


    use chrono::*;
    use typenum::U3;

    #[test]
    fn bench_bb_alg() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]);
        println!("BB Alg time = {:?}", 
            Duration::span(|| {
                r.bb_algorithm();
            }));

        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
                Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
            ];

        let ideal = Ideal::new(f_vec);

        println!("BB Alg time = {:?}", 
            Duration::span(|| {
                ideal.bb_algorithm();
            }));
    }

    // #[test]
    // fn bench_bb_alg_impr() {
    //     let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
    //     let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
    //     let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
    //     let r = &Ideal::new(vec![a, b]);
    //     println!("BB Alg time = {:?}", 
    //         Duration::span(|| {
    //             r.bb_algorithm_impr();
    //         }));

    //     let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
    //     let f_vec = vec![
    //             Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
    //             Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
    //             Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
    //         ];

    //     let ideal = Ideal::new(f_vec);

    //     println!("BB Alg time = {:?}", 
    //         Duration::span(|| {
    //             ideal.bb_algorithm_impr();
    //         }));
    // }

    /// Tests the "is_groebner_basis()" function
    #[test]
    fn is_gb_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let x = Poly::from_str(&ring, "1.0x^1").unwrap();
        let y = Poly::from_str(&ring, "1.0y^1").unwrap();
        let x_y_ideal = Ideal::new(vec![x, y]);
        assert!(x_y_ideal.is_groebner_basis());

        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let not_gb = Ideal::new(vec![a, b]);
        assert!(!not_gb.is_groebner_basis());
    }
}
