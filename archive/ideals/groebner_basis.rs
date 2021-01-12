use crate::ideals::*;

// use crate::algebras::EuclideanDomain;

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
use num_traits::Zero;

impl<'a, P: FPolyRing> Ideal<'a, P> {
    /// Standard implementation of Buchberger's Criterion to check if a
    /// basis is a Groebner Basis
    pub fn is_groebner_basis(&self) -> bool {
        let n = self.gens.len();
        iproduct!(0..n, 0..n)
            .filter(|(i, j)| i <= j)
            .map(|(i, j)| self.gens[i].s_poly(&self.gens[j]).unwrap().reduce(&self))
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
                let r = gb.gens[i].s_poly(&gb.gens[j]).unwrap().reduce(&gb);
                if !r.is_zero() {
                    // Add the remainder to the ideal and update the pairs that we
                    // need to check
                    if let Some(_) = gb.add(r) {
                        n += 1;
                        new_pairs.append(&mut (0..n - 1).map(|k| (k, n - 1)).collect());
                    }
                }
            }
            pairs = new_pairs;
        }
        gb
    }

    //     /// Implementation of the Improved Buchberger's algorithm found at the end of
    //     /// "Ideals, Varieties, and Algorithms" by Cox, Little, and O'Shea 4th edition.
    //     /// Not fully implemented, the trick from Proposition 3 hasn't yet been implemented as it is
    //     /// non-trivial
    // pub fn bb_algorithm_impr(&self) -> Self {
    //     let mut gb = self.clone();
    //     let mut n = gb.gens.len();
    //     let mut unchecked_pairs = unord_pairs_int(n);

    //     // Iterates through all the pairs until there are no more pairs to check
    //     while !unchecked_pairs.is_empty() {
    //         let mut new_pairs = Vec::new();

    //         for (i, j) in unchecked_pairs {
    //             if gb.gens[i].lt().gcd(&gb.gens[j].lt()) != Term::one() {
    //                 let r = s_poly(&gb.gens[i], &gb.gens[j]).reduce(&gb);
    //                 if !r.is_zero() {
    //                     // Add the remainder to the ideal and update the pairs that we
    //                     // need to check
    //                     gb.add(r);
    //                     n += 1;
    //                     new_pairs.append(&mut (0..n - 1).map(|k| (k, n)).collect());
    //                 }
    //             }
    //         }
    //         unchecked_pairs = new_pairs;
    //     }
    //     gb
    // }
}

fn unord_pairs_int(n: usize) -> Vec<(usize, usize)> {
    iproduct!((0..n), (0..n)).filter(|(i, j)| i >= j).collect()
}

impl<'a, P: FPolyRing> Poly<'a, P> {
    /// The remainder upon division by the polynomials in "g".
    /// Note: If "g" is not a Groebner basis then this can produce unpredictable
    /// results
    fn reduce(&self, g: &Ideal<'a, P>) -> Self {
        self.divpolys(&g.gens).1
    }

    fn s_poly(&self, rhs: &Poly<'a, P>) -> Option<Poly<'a, P>> {
        if self.is_zero() || rhs.is_zero() {
            None
        } else {
            let lcm = self.lt().lcm(&rhs.lt());
            let a_newlead = lcm.euclid_div(&self.lt()).unwrap().0;
            let b_newlead = lcm.euclid_div(&rhs.lt()).unwrap().0;
            Some((self * a_newlead) - (rhs * b_newlead))
        }
    }

    /// Divides self by the vector of polynomial in the order they were given in
    /// This is an implementation of the algorithm given in CLO
    pub fn divpolys(&self, divisors: &Vec<Self>) -> (Vec<Self>, Self) {
        // Since we need to mutate self. There is a potential optimisation to be had by
        // getting rid of this clone
        let mut p = self.clone();

        // Remainder and quotients accumulators
        let mut r = Poly::zero();
        let mut q = vec![Poly::zero(); divisors.len()];

        'outer: while !p.is_zero() {
            for i in 0..divisors.len() {
                // Checks if the lead term of p is divisible by one of the values
                let quo = p.lt().euclid_div(&divisors[i].lt()).unwrap().0;
                if !quo.is_zero() {
                    p = p - (&divisors[i] * quo.clone());
                    q[i] = &q[i] + &Poly::from_terms(vec![quo], q[i].ring);
                    continue 'outer;
                }
            }
            // Executes if no division occurred
            r = r + Poly::from_terms(vec![p.terms.pop().unwrap()], self.ring);
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
        println!("1");

        let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);

        // Multivariate
        // let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        // let a = Poly::from_str(&ring, "3.0x^2y^6 + 5.0x^98y^2").unwrap();
        // let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        // println!("a = {}", a);
        // println!("b = {}", b);
        // let (q, r) = a.divpolys(&vec![b]);
        // println!("q = {:?}", q);
        // println!("r = {}", r);
    }

    /// Tests the normal Buchberger's Algorithm
    #[test]
    fn bb_alg_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]).bb_algorithm();
        assert!(r.is_groebner_basis());
    }

    use chrono::*;
    use typenum::U3;

    #[test]
    fn bench_bb_alg() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]);
        println!(
            "BB Alg time = {:?}",
            Duration::span(|| {
                r.bb_algorithm();
            })
        );

        let ring = PRDomain::<RR, GLex<MultiIndex<U3>>>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
            Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
            Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
            Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
        ];

        let ideal = Ideal::new(f_vec);

        println!(
            "BB Alg time = {:?}",
            Duration::span(|| {
                ideal.bb_algorithm();
            })
        );
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
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let x = Poly::from_str(&ring, "1.0x^1").unwrap();
        let y = Poly::from_str(&ring, "1.0y^1").unwrap();
        let x_y_ideal = Ideal::new(vec![x, y]);
        assert!(x_y_ideal.is_groebner_basis());
        println!("waypoint");

        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let not_gb = Ideal::new(vec![a, b]);
        assert!(!not_gb.is_groebner_basis());
    }
}
